%% Over-approximation of the sensitivity bounds using second-order sensitivities
% Steps:
%   1. Over-approximate reachable tube from t_init to t_final of 
%       first-order sensitivity matrices.
%   2. Over-approximate the reachable set at time t_final of the 
%       second-order sensitivity matrices using the reachable tube
%       of first-order sensitivity matrices.
%   3. Use quasi-Monte Carlo methods similar to OA_6_DT_CT_Quasi_Monte_Carlo.m
%       to grid-sample the initial states and inputs at t_init, evolve
%       these in time to t_final, and compute first-order sensitivity matrices.
%   4. Create an interval approximation of the first-order sensitivity
%       matrices from step 3, and expand to an over-approximation
%       using an error-bound generated from the second-order sensitivity matrices.
%   5. Return this first-order sensitivity matrices interval over-approximation.

% This method guarantees an over-approximation of the sensitivity
% matrices intervals, unlike the sampling and falsification method.

% This method can be used to generate arbitrarily tight intervals
% by increasing the number of samples used in samples_per_state_dimension
% and samples_per_input dimension, since as the number of samples increases,
% the dispersion decreases and goes to 0, making the error bound go to 0.
% Increasing the number of samples comes at computational cost, so
% the number of samples per dimension can be used to tune between
% computation time and interval tightness.
% This method may produce significantly less conservative estimates
% than the interval arithmetics method.

% Source paper:
% P.-J. Meyer and M. Arcak, "Interval Reachability Analysis using 
% Second-Order Sensitivity". Proceedings of the 21st IFAC World Congress, 
% pp. 1851-1856, 2020. DOI: 10.1016/j.ifacol.2020.12.2344

% List of inputs:
%   t_init: initial time
%   t_final: time for which the sensitivity matrices interval 
%       over-approximation is computed.
%   [x_low, x_up]: interval of initial states (at time t_init).
%   [p_low, p_up]: interval of allowed input values.
%   samples_per_state_dimension: vector of same length as state variable. 
%       The i'th value of this vector is the number of 
%       points to sample along dimension i of the initial state interval.
%   samples_per_input_dimension: vector of same length as input variable.
%       The j'th value of this vector is the number of points
%       to sample along dimension j of the interval over input values.
%   sensitivity_equation_solver: 
%       0 Euler method on the sensitivity definition
%       1 Euler method on Jacobian + ODE solver on sensitivity ODE
%       2 ODE solver on sensitivity ODE using explicit Jacobian function (requires Input_files/UP_Jacobian_Function.m)
%       See Sensitivity_solver.m solver_choice for details.
%   [J_*_low, J_*_up]: Jacobian interval bounds.
%   [J_**_low, J_**_up]: Second-order Jacobian interval bounds.

% List of outputs:
%   [S_x_low, S_x_up]: Sensitivity matrix interval over-approximation
%       with respect to initial state.
%   [S_p_low, S_p_up]: Sensitivity matrix interval over-approximation
%       with respect to input.

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [S_x_low, S_x_up, S_p_low, S_p_up] = Sensitivity_from_second_order(...
    t_init, t_final, x_low, x_up, p_low, p_up,...
    samples_per_state_dimension, samples_per_input_dimension,...
    sensitivity_equation_solver,...
    J_x_low, J_x_up, J_p_low, J_p_up,...
    J_xx_low, J_xx_up, J_xp_low, J_xp_up,...
    J_pp_low, J_pp_up, J_px_low, J_px_up)

import interval_arithmetics.*;
% Alias some function names for conciseness.
rstp = @Matrix_right_semitensor_product;
rstp_ia = @Interval_right_semitensor_product;
lstp_ia = @Interval_left_semitensor_product;

log = logger.get_logger();

assert(isscalar(t_init), 't_init must be a scalar');
assert(isscalar(t_final), 't_final must be a scalar');

% Verify size and ordering of initial state interval.
assert(all(size(x_low) == size(x_up)),...
    'size of x_low must match size of x_up');
assert(all(x_up >= x_low),...
    'x_up must be greater than or equal to x_low, component-wise');

% Verify size and ordering of input interval.
assert(all(size(p_low) == size(p_up)),...
    'size of p_low must match size of p_up');
assert(all(p_up >= p_low),...
    'p_up must be greater than or equal to p_low, component-wise');

% Verify number of points per state dimension.
assert(length(samples_per_state_dimension) == length(x_low),...
    'size of samples_per_state_dimension must match size of state x_low');
assert(all(samples_per_state_dimension >= 1),...
    'The number of samples per state dimension must be a positive integer');
assert(all(mod(samples_per_state_dimension, 1) == 0),...
    'The number of samples per state dimension must be a positive integer');

% Verify number of points per input dimension.
assert(length(samples_per_input_dimension) == length(p_low),...
    'size of samples_per_input_dimension must match size of input p_low');
assert(all(samples_per_input_dimension >= 1),...
    'The number of samples per input dimension must be a positive integer');
assert(all(mod(samples_per_input_dimension, 1) == 0),...
    'The number of samples per input dimension must be a positive integer');

assert(isscalar(sensitivity_equation_solver),...
    'sensitivity_equation_solver must be a scalar');

n_x = numel(x_low);
n_p = numel(p_low);

%% Initialization
% Ignore state or input sensitivity if the corresponding interval is a singleton
needed_sensitivies = 0;         % Default: need all 6 sensitivities
if isequal(x_low,x_up)
    needed_sensitivies = 2;     % Only compute S_p and S_pp
    S_x_low = zeros(n_x);
    S_x_up = zeros(n_x);
    S_px_low = zeros(n_x,n_x*n_p);
    S_px_up = zeros(n_x,n_x*n_p);
elseif isequal(p_low,p_up)
    needed_sensitivies = 1;     % Only compute S_x and S_xx
    S_p_low = zeros(n_x,n_p);
    S_p_up = zeros(n_x,n_p);
    S_xp_low = zeros(n_x,n_x*n_p);
    S_xp_up = zeros(n_x,n_x*n_p);
end

tic_2ndSensi = tic;

%% Over-approximate the reachable TUBE for first-order sensitivities
% using interval arithmetics

% First-order sensitivity with respect to initial states
if needed_sensitivies ~= 2
    [S_x_RT_low,S_x_RT_up] = Affine_interval_system_reachability(t_init,t_final,...
        eye(n_x),eye(n_x),J_x_low,J_x_up,zeros(n_x),zeros(n_x),1);
end

% First-order sensitivity with respect to inputs
if needed_sensitivies ~= 1
    [S_p_RT_low,S_p_RT_up] = Affine_interval_system_reachability(t_init,t_final,...
        zeros(n_x,n_p),zeros(n_x,n_p),J_x_low,J_x_up,J_p_low,J_p_up,1);
end

%% Over-approximate the reachable SET for second-order sensitivities
% using interval arithmetics

% Sensitivity with respect to initial states twice
if needed_sensitivies ~= 2
    [temp1_LB,temp1_UB] = rstp_ia(J_xx_low,J_xx_up,S_x_RT_low,S_x_RT_up);                           % J_xx*S_x with right STP
    [temp2_LB,temp2_UB] = lstp_ia(temp1_LB,temp1_UB,S_x_RT_low,S_x_RT_up);                          % (J_xx*S_x)*S_x with left STP
    [S_xx_low,S_xx_up] = Affine_interval_system_reachability(t_init,t_final,...
        zeros(n_x,n_x^2),zeros(n_x,n_x^2),J_x_low,J_x_up,temp2_LB,temp2_UB,0);
end

% Sensitivity with respect to initial states then inputs
if needed_sensitivies == 0
    [temp1_LB,temp1_UB] = rstp_ia(J_xx_low,J_xx_up,S_p_RT_low,S_p_RT_up);                   % J_xx*S_p with right STP
    [temp2_LB,temp2_UB] = Interval_matrix_sum(temp1_LB,temp1_UB,J_xp_low,J_xp_up);          % J_xx*S_p + J_xp
    [temp3_LB,temp3_UB] = lstp_ia(temp2_LB,temp2_UB,S_x_RT_low,S_x_RT_up);                  % (J_xx*S_p + J_xp)*S_x with left STP
    [S_xp_low,S_xp_up] = Affine_interval_system_reachability(t_init,t_final,...
        zeros(n_x,n_x*n_p),zeros(n_x,n_x*n_p),J_x_low,J_x_up,temp3_LB,temp3_UB,0);
end

% Sensitivity with respect to inputs twice
if needed_sensitivies ~= 1
    [temp1_LB,temp1_UB] = rstp_ia(J_xx_low,J_xx_up,S_p_RT_low,S_p_RT_up);                   % J_xx*S_p with right STP
    [temp2_LB,temp2_UB] = Interval_matrix_sum(temp1_LB,temp1_UB,J_xp_low,J_xp_up);          % J_xx*S_p + J_xp
    [temp3_LB,temp3_UB] = lstp_ia(temp2_LB,temp2_UB,S_p_RT_low,S_p_RT_up);                  % (J_xx*S_p + J_xp)*S_p with left STP
    [temp4_LB,temp4_UB] = rstp_ia(J_px_low,J_px_up,S_p_RT_low,S_p_RT_up);                   % J_px*S_p with right STP
    [temp5_LB,temp5_UB] = Interval_matrix_sum(temp4_LB,temp4_UB,J_pp_low,J_pp_up);          % J_px*S_p + J_pp
    [temp6_LB,temp6_UB] = Interval_matrix_sum(temp3_LB,temp3_UB,temp5_LB,temp5_UB);         % (J_xx*S_p + J_xp)*S_p + J_px*S_p + J_pp
    [S_pp_low,S_pp_up] = Affine_interval_system_reachability(t_init,t_final,...
        zeros(n_x,n_p^2),zeros(n_x,n_p^2),J_x_low,J_x_up,temp6_LB,temp6_UB,0);
end

% Sensitivity with respect to inputs then initial states
if needed_sensitivies == 0
    [temp1_LB,temp1_UB] = rstp_ia(J_xx_low,J_xx_up,S_x_RT_low,S_x_RT_up);                           % J_xx*S_x with right STP
    [temp2_LB,temp2_UB] = lstp_ia(temp1_LB,temp1_UB,S_p_RT_low,S_p_RT_up);                          % (J_xx*S_x)*S_p with left STP
    [temp3_LB,temp3_UB] = rstp_ia(J_px_low,J_px_up,S_x_RT_low,S_x_RT_up);                           % J_px*S_x with right STP
    [temp4_LB,temp4_UB] = Interval_matrix_sum(temp2_LB,temp2_UB,temp3_LB,temp3_UB);                 % (J_xx*S_x)*S_p + J_px*S_x
    [S_px_low,S_px_up] = Affine_interval_system_reachability(t_init,t_final,...
        zeros(n_x,n_x*n_p),zeros(n_x,n_x*n_p),J_x_low,J_x_up,temp4_LB,temp4_UB,0);
end

%% Over-approximate the reachable SET for first-order sensitivities
% by evaluating the first-order sensitivities for a few samples of the
% initial states and inputs, then using the bounds on the second-order
% sensitivities to obtained guaranteed bounds on the first-order ones.

% Number of samples per dimension
if needed_sensitivies == 2
    samples_per_state_dimension = ones(n_x,1);
end
if needed_sensitivies == 1
    samples_per_input_dimension = ones(n_p,1);
end
sensi_sampling_total = prod(samples_per_state_dimension)*prod(samples_per_input_dimension);

% Convert sample indices to subscripts on each dimension
x0_subscripts = cell(n_x,1);
[x0_subscripts{:}] = ind2sub(samples_per_state_dimension',1:prod(samples_per_state_dimension));
x0_subscripts = cell2mat(x0_subscripts);
p_subscripts = cell(n_p,1);
[p_subscripts{:}] = ind2sub(samples_per_input_dimension',1:prod(samples_per_input_dimension));
p_subscripts = cell2mat(p_subscripts);

% Compute componentwise dispersion for a uniform grid
%   the step between two samples is twice the dispersion
dispersion_x = (x_up-x_low)./(samples_per_state_dimension*2);
dispersion_p = (p_up-p_low)./(samples_per_input_dimension*2);

% Evaluate the first-order sensitivities
S_x_samples = NaN(n_x^2,sensi_sampling_total);
S_p_samples = NaN(n_x*n_p,sensi_sampling_total);
for x0_ind = 1:prod(samples_per_state_dimension)
    x0_samples = x_low + dispersion_x + 2*dispersion_x.*(x0_subscripts(:,x0_ind)-1);
    for p_ind = 1:prod(samples_per_input_dimension)
        p_samples = p_low + dispersion_p + 2*dispersion_p.*(p_subscripts(:,p_ind)-1);
        i = sub2ind([prod(samples_per_state_dimension) prod(samples_per_input_dimension)],x0_ind,p_ind);

        % Compute corresponding sensitivity values
        [S_x,S_p] = Sensitivity_solver(t_init,t_final,x0_samples,p_samples,...
            needed_sensitivies,sensitivity_equation_solver);
        S_x_samples(:,i) = S_x(:);
        S_p_samples(:,i) = S_p(:);
    end
end

% Convert to scalar dispersion and then form vectors
dispersion_x = ones(n_x,1)*norm(dispersion_x,Inf);
dispersion_p = ones(n_p,1)*norm(dispersion_p,Inf);

% Get the extrema of the sensitivitiy evaluations
S_x_sample_min = reshape(min(S_x_samples,[],2),n_x,n_x);
S_x_sample_max = reshape(max(S_x_samples,[],2),n_x,n_x);
S_p_sample_min = reshape(min(S_p_samples,[],2),n_x,n_p);
S_p_sample_max = reshape(max(S_p_samples,[],2),n_x,n_p);

% First-order sensitivity with respect to initial states
if needed_sensitivies ~= 2
    S_x_low = S_x_sample_min - rstp(max(abs(S_xx_low),abs(S_xx_up)),dispersion_x) - rstp(max(abs(S_xp_low),abs(S_xp_up)),dispersion_p);
    S_x_up = S_x_sample_max + rstp(max(abs(S_xx_low),abs(S_xx_up)),dispersion_x) + rstp(max(abs(S_xp_low),abs(S_xp_up)),dispersion_p);
end

% First-order sensitivity with respect to inputs
if needed_sensitivies ~= 1
    S_p_low = S_p_sample_min - rstp(max(abs(S_px_low),abs(S_px_up)),dispersion_x) - rstp(max(abs(S_pp_low),abs(S_pp_up)),dispersion_p);
    S_p_up = S_p_sample_max + rstp(max(abs(S_px_low),abs(S_px_up)),dispersion_x) + rstp(max(abs(S_pp_low),abs(S_pp_up)),dispersion_p);
end

toc_2ndSensi = toc(tic_2ndSensi);
log.runtime('Computing sensitivity bounds took %f seconds\n', toc_2ndSensi);
end % end Sensitivity_from_second_order
