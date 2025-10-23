%% Compute the reachable set or tube of an affine interval system
% dz/dt = A*z + B, with matrices A and B bounded in intervals

% Source paper 1 (interval arithmetics for time-varying linear systems):
% M. Althoff, O. Stursberg and M. Buss, "Reachability analysis of linear 
% systems with uncertain parameters and inputs". 46th IEEE Conference on 
% Decision and Control, pp. 726-732, 2007. DOI: 10.1109/CDC.2007.4434084

% Source paper 2 (definition of basic operations for interval arithmetics):
% L. Jaulin, M. Kieffer, O. Didrit and E. Walter, "Applied interval 
% analysis: with examples in parameter and state estimation, robust control
% and robotics". Springer Science & Business Media, v. 1, 2001.

% List of inputs
%   t_init: initial time
%   t_final: time at which the reachable set is approximated
%   [z0_LB,z0_UB]: interval of initial states, with size(z) = [n1,n2]
%   [A_LB,A_UB]: bounds of the Jacobian matrix A, with size(A) = [n1,n1]
%   [B_LB,B_UB]: bounds of the input matrix B, with size(B) = [n1,n2]
%   bool_reachable_tube: choice of what to over-approximate:
%       0 for the reachable set at t_final
%       1 for the reachable tube over [t_init,t_final]
%   taylor_order: (optional) integer to impose a minimum order of the Taylor series
%   scaling_factor: (optional) scaling factor used in Scaling and Squaring 
%       method for computation of exponential series

% List of outputs
%   [z_LB,z_UB]: interval over-approximation of the reachable set or tube

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function [z_LB,z_UB] = Affine_interval_system_reachability(t_init, t_final,...
    z0_LB, z0_UB, A_LB, A_UB, B_LB, B_UB, bool_reachable_tube, taylor_order, scaling_factor)

% Dummy state and input dimensions to allow Solver_parameters.m to run properly
n_x = 1;
n_p = 1;

import interval_arithmetics.*;

[n1,n2] = size(z0_LB);
T = t_final-t_init;     % Width of the time step

% Check consistency of user-provided inputs
assert(t_init<t_final,'Inconsistent time interval: need t_init<t_final')
assert(isequal(size(z0_LB),size(z0_UB)),'Inconsistent size of the interval of initial states: size(z0_LB)~=size(z0_UB)')
assert(isequal(size(A_LB),[n1,n1]),'Inconsistent size of the Jacobian interval matrix: need A_LB and A_UB to be square of dimension size(z0_LB,1)')
assert(isequal(size(A_LB),size(A_UB)),'Inconsistent size of the Jacobian interval matrix: size(A_LB)~=size(A_UB)')
assert(isequal(size(B_LB),size(z0_UB)),'Inconsistent size of the input interval matrix: need B_LB and B_UB to have the same size as z0_LB')
assert(isequal(size(B_LB),size(B_UB)),'Inconsistent size of the input interval matrix: size(B_LB)~=size(B_UB)')
assert(all(z0_LB(:)<=z0_UB(:)),'Inconsistent bounds of the initial state: need z0_LB <= z0_UB')
assert(all(A_LB(:)<=A_UB(:)),'Inconsistent bounds of the Jacobian matrix: need A_LB <= A_UB')
assert(all(B_LB(:)<=B_UB(:)),'Inconsistent bounds of the input matrix: need B_LB <= B_UB')

% Extract from Solver_parameters.m the parameters for Taylor series truncation
run('Solver_parameters.m')
% Only read 'scaling_factor' and 'taylor_order' from Solver_parameters.m if not provided as an input to this function.
if nargin < 11
    scaling_factor = parameters.scaling_factor;     % Factor by which to scale down matrix A*T to improve convergence of homogeneous solution. 
end
if nargin < 10 
    taylor_order = parameters.taylor_order;         % Initial order for the truncation of the Taylor series (may be interatively increased if the results are not satisfactory)
end
taylor_tolerance = parameters.taylor_tolerance;     % Stopping condition: max relative error allowed between two iterations
taylor_increment = parameters.taylor_increment;     % Taylor order increase at each iteration

% Increase Taylor order if it is too low for the computations to hold
taylor_order = max(3, max(taylor_order, floor(Interval_matrix_norm(A_LB, A_UB)*T-2) + 1));

%% Iterative computation of the reachable set with increasing Taylor orders until convergence
% To skip this (do a single iteration), set: 'parameters.taylor_tolerance = NaN' in Solver_parameters.m 

% Initial computation of the reachable set
[z_LB, z_UB] = Reachable_set_OA(T, taylor_order, scaling_factor, z0_LB, z0_UB, A_LB, A_UB, B_LB, B_UB, bool_reachable_tube);

% Initialization of the internal variables for the loop update
z_LB_prev = zeros(n1,n2);
z_UB_prev = zeros(n1,n2);

% Loop stops when  interval [z_LB,z_UB] has not changed more than tolerance
% or if taylor_tolerance is undefined (NaN)
while ~isnan(taylor_tolerance) && ...
      (~all(all((abs(z_LB_prev-z_LB) <= abs(z_LB)*taylor_tolerance) & ...
                (abs(z_UB_prev-z_UB) <= abs(z_UB)*taylor_tolerance))))
    % Update previous reachable set/tube
    z_LB_prev = z_LB;
    z_UB_prev = z_UB;
    
    % Increase Taylor order
    taylor_order = taylor_order + taylor_increment;
    
    % Compute new reachable set/tube
    [z_LB, z_UB] = Reachable_set_OA(T, taylor_order, scaling_factor, z0_LB, z0_UB, A_LB, A_UB, B_LB, B_UB, bool_reachable_tube);
end
end % End Affine_interval_system_reachability.


%% Computation of the reachable set or tube of the affine interval system
% Using truncated Taylor series and an over-approximation of its remainder

% This function computes, for system dz/dt=Az+B:
% - an over-approximation of the homogeneous solution after time T
% - an over-approximation of the particulate solution after time T
% - an over-approximation of the homogeneous solution over time range [0,T]
% - an over-approximation of the reachable set after time T by adding both
% homogeneous and particulate solutions after T
% - or an over-approximation of the reachable tube over time range [0,T] by
% adding the homogeneous solution over time range [0,T] with the
% particulate solution at time T.

% List of inputs
%   T: time step
%   taylor_order: truncation order of the Taylor series
%   scaling_factor: scaling factor used in Scaling and Squaring method for
%       computation of exponential series.
%   [z0_LB,z0_UB]: interval of initial states
%   [A_LB,A_UB]: bounds of the Jacobian matrix A
%   [B_LB,B_UB]: bounds of the input matrix B
%   bool_reachable_tube: choice of what to over-approximate:
%       0 for the reachable set at t_final
%       1 for the reachable tube over [t_init,t_final]

% List of outputs
%   [z_LB,z_UB]: interval over-approximation of the reachable set or tube
function [z_LB, z_UB] = Reachable_set_OA(T, taylor_order, scaling_factor, z0_LB, z0_UB, A_LB, A_UB, B_LB, B_UB, bool_reachable_tube)

import interval_arithmetics.*;

[n1,n2] = size(z0_LB);

% Particular case if A is a singleton equal to the 0 matrix
if isequal(A_LB,A_UB) && ~any(A_LB(:))
    % Reachable set at the final time
    [z_LB,z_UB] = Interval_matrix_sum(z0_LB,z0_UB,B_LB*T,B_UB*T);
    
    % Reachable tube if necessary
    if bool_reachable_tube
        % Interval hull with the set of initial states
        z_LB = min(z0_LB,z_LB);
        z_UB = max(z0_UB,z_UB);
    end
    return
end

% Initialization of a variables that are reused a lot
X_LB = A_LB*T;
X_UB = A_UB*T;
scalar_LB = @(i) i^(i/(1-i)) - i^(1/(1-i));

% Remainder of the Taylor series after truncation at taylor_order
[res_LB,res_UB] = Truncation_remainder(X_LB,X_UB,taylor_order);

%% Homogeneous solution (and sum component of the hull enlargement term)
if isequal(z0_LB,z0_UB) && ~any(z0_LB(:))
    % Skipped if the initial state interval is a singleton equal to 0
    hom_LB = zeros(n1,n2);
    hom_UB = zeros(n1,n2);
    HE_LB = zeros(n1,n2);
    HE_UB = zeros(n1,n2);
else
    % Over-approximation of the interval matrix exponential exp(A*T)
    [expAT_LB,expAT_UB] = Interval_matrix_exponential(X_LB,X_UB,taylor_order,scaling_factor);
    
    % Reachable set component from the homogeneous solution
    [hom_LB,hom_UB] = Interval_matrix_product(expAT_LB,expAT_UB,z0_LB,z0_UB);

    % Additional computations if the reachable TUBE is also needed
    if bool_reachable_tube
        % Compute truncated sum of [scalar_LB(i),0]*(A*T)^i/i!
        % from i=2 to i=taylor_order, using factored form
        
        % Initialization of the hull enlargement term due to initial states
        % [scalar_LB(taylor_order),0]*A*T
        [HE_LB,HE_UB] = Interval_matrix_product_with_scalar(scalar_LB(taylor_order),0,X_LB,X_UB);
        
        % Loop of: HE = HE*A*T/i + [scalar_LB(i-1),0]*A*T
        for i = taylor_order:-1:2
            [HE_LB,HE_UB] = Interval_matrix_product(HE_LB,HE_UB,X_LB/i,X_UB/i);
            [temp_LB,temp_UB] = Interval_matrix_product_with_scalar(scalar_LB(i-1),0,X_LB,X_UB);
            [HE_LB,HE_UB] = Interval_matrix_sum(HE_LB,HE_UB,temp_LB,temp_UB);
        end

        % Addition of the remainder from the Taylor series truncation
        [HE_LB,HE_UB] = Interval_matrix_sum(HE_LB,HE_UB,res_LB,res_UB);

        % Final hull enlargement term due to initial states
        [HE_LB,HE_UB] = Interval_matrix_product(HE_LB,HE_UB,z0_LB,z0_UB);
    end
end

%% Particulate solution
if isequal(B_LB,B_UB) && ~any(B_LB(:))
    % Skipped if the input interval is a singleton equal to 0
    par_LB = zeros(n1,n2);
    par_UB = zeros(n1,n2);

    par_center_LB = zeros(n1,n2);
    par_center_UB = zeros(n1,n2);
    HE_input_LB = zeros(n1,n2);
    HE_input_UB = zeros(n1,n2);
else
    % Check if 0 belongs to the input set [B_LB,B_UB] (for reachable tube)
    if bool_reachable_tube && ~(all(B_LB(:) <= zeros(n1*n2,1)) && all(B_UB(:) >= zeros(n1*n2,1)))
        % If not, redefine [B_LB,B_UB] := [B_LB,B_UB] - B_center
        % so that 0 belongs to the new [B_LB,B_UB]
        B_center = (B_LB+B_UB)/2;
        [B_LB,B_UB] = Interval_matrix_sum(B_LB,B_UB,-B_center,-B_center);

        bool_zero_in_B = 0;     % This also implies bool_reachable_tube is true
    else
        bool_zero_in_B = 1;     % This might mean that bool_reachable_tube is false
        % Set unnecessary variables to 0
        par_center_LB = zeros(n1,n2);
        par_center_UB = zeros(n1,n2);
        HE_input_LB = zeros(n1,n2);
        HE_input_UB = zeros(n1,n2);
    end

    % The particulate solution (reachable set component with respect to the
    % input set) relies on the sum of (A*T)^i/(i+1)!
    % that we compute below in the factored form
    par_LB = eye(n1);
    par_UB = eye(n1);
    % Loop of: par = par*A*T/(i+1) + eye(n1)
    for i = taylor_order:-1:1
        [par_LB,par_UB] = Interval_matrix_product(par_LB,par_UB,X_LB/(i+1),X_UB/(i+1));
        [par_LB,par_UB] = Interval_matrix_sum(par_LB,par_UB,eye(n1),eye(n1));
    end
    
    % Addition of the remainder from the Taylor series truncation
    [par_LB,par_UB] = Interval_matrix_sum(par_LB,par_UB,res_LB,res_UB);
    
    % Reachable set component from the center of the input set
    %   (only needed if the input set does not contain 0)
    if ~bool_zero_in_B
        [par_center_LB,par_center_UB] = Interval_matrix_product(par_LB,par_UB,B_center*T,B_center*T);
    end
    
    % Reachable set component from the main or shifted input set
    [par_LB,par_UB] = Interval_matrix_product(par_LB,par_UB,B_LB*T,B_UB*T);

    % Additional computations of a hull enlargement if the reachable TUBE 
    % is also needed and if the input set does not contain 0
    if ~bool_zero_in_B
        % Compute truncated sum of [scalar_LB(i),0]*(A*T)^(i-1)/i!
        % from i=2 to i=taylor_order, using factored form
        
        % Initialization of the hull enlargement term due to the shifted
        % input center B_center: [scalar_LB(taylor_order),0]*A*T
        [HE_input_LB,HE_input_UB] = Interval_matrix_product_with_scalar(scalar_LB(taylor_order),0,X_LB,X_UB);
        
        % Loop of: HE_input = HE_input*A*T/i + [scalar_LB(i-1),0]*A*T
        for i = taylor_order:-1:2
            [HE_input_LB,HE_input_UB] = Interval_matrix_product(HE_input_LB,HE_input_UB,X_LB/i,X_UB/i);
            [temp_LB,temp_UB] = Interval_matrix_product_with_scalar(scalar_LB(i-1),0,X_LB,X_UB);
            [HE_input_LB,HE_input_UB] = Interval_matrix_sum(HE_input_LB,HE_input_UB,temp_LB,temp_UB);
        end
        
        % Adapt result to obtain the sum of [scalar_LB(i),0]*(A*T)^(i-1)/i!
        [HE_input_LB,HE_input_UB] = Interval_matrix_product_with_scalar(T/2,T/2,HE_input_LB,HE_input_UB);
        
        % Add error term for the input-related hull enlargement
        inf_norm = Interval_matrix_norm(A_LB,A_UB);
        [HE_input_LB,HE_input_UB] = Interval_matrix_sum(HE_input_LB,HE_input_UB,res_LB/inf_norm,res_UB/inf_norm);
        
        % Final input-related hull enlargement
        [HE_input_LB,HE_input_UB] = Interval_matrix_product(HE_input_LB,HE_input_UB,B_center,B_center);
    end
end

%% Resulting reachable set or reachable tube
if ~bool_reachable_tube
    % Final computations for the reachable set
    [z_LB,z_UB] = Interval_matrix_sum(hom_LB,hom_UB,par_LB,par_UB);
else
    % Final computations for the reachable tube
    
    % Homogeneous solution + influence of the input center
    [z_LB,z_UB] = Interval_matrix_sum(hom_LB,hom_UB,par_center_LB,par_center_UB);

    % Interval hull with the set of initial states
    z_LB = min(z0_LB,z_LB);
    z_UB = max(z0_UB,z_UB);

    % Enlargement of the hull related to the homogeneous solution
    [z_LB,z_UB] = Interval_matrix_sum(z_LB,z_UB,HE_LB,HE_UB);

    % Enlargement of the hull related to the input center
    [z_LB,z_UB] = Interval_matrix_sum(z_LB,z_UB,HE_input_LB,HE_input_UB);
    
    % Final reachable tube (after adding particulate solution)
    [z_LB,z_UB] = Interval_matrix_sum(z_LB,z_UB,par_LB,par_UB);
end
end % End Reachable_set_OA.


%% Over-approximation of the remainder of the truncated Taylor series
% Inputs:
%   [A_low, A_up]: interval matrices A.
%   taylor_order: truncation order of the Taylor series.
% Outputs:
%   [E_low, E_up]: bounds of the remainder of the truncated Taylor series
function [E_low,E_up] = Truncation_remainder(A_low, A_up, taylor_order)
    import interval_arithmetics.Interval_matrix_norm;

    [n_x,~] = size(A_low);
    A_norm = Interval_matrix_norm(A_low,A_up);
    E_up = (A_norm^(taylor_order + 1)) / (factorial(taylor_order + 1) * (1 - A_norm / (taylor_order + 2)));
    E_up = ones(n_x)*E_up;
    E_low = -E_up;
end % End Truncation_remainder.
