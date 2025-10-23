%% Compute interval matrices bounding the possible values of each sensitivity matrix
% The user can either ask for a specific method, or let the function pick
% the best method based on which user-provided files are available in the
% folder ./Input_files
% For the sampling and falsification submethod and second-order sensitivity
% submethod, the user can also specify which solver to use to evaluate 
% the sensitivity values.

% Source paper:
% P.-J. Meyer, A. Devonport and M. Arcak, "Sampled-Data Mixed Monotonicity" 
% In: "Interval Reachability Analysis: Bounding Trajectories of Uncertain 
% Systems with Boxes for Control and Verification". Springer, 2021
% DOI: 10.1007/978-3-030-65110-7_5

% List of inputs
%   t_init: initial time.
%   t_final: time at which the reachable set is approximated.
%   [x_low,x_up]: interval of initial states (at time t_init).
%   [p_low,p_up]: interval of allowed input values.

% Optional inputs (read from Solver_parameters.m if not provided in the call)
%   sensitivity_bounds_method: integer to pick the submethod to obtain the sensitivity bounds
%       1: Read sensitivity bounds from user-provided files 
%           Requires either Input_files/UP_Sensitivity_Bounds.m 
%           or Input_files/UP_Sensitivity_Signs.m
%       2: Interval Arithmetics 
%           Requires Input_files/UP_Jacobian_Function.m
%       3: 2nd order Sensitivity 
%           Requires Input_files/UP_Jacobian_Bounds.m 
%           and Input_files/UP_Jacobian_2ndOrder_Bounds.m
%       0: Sampling and falsification 
%           No additional requirement if sensitivity_equation_solver~=2,
%           else requires Input_files/UP_Jacobian_Function.m
%       NaN: let this function select the best available method based on 
%           provided files in './Input_files' folder
%           (Default if no file is provided: 0) 
%   sensitivity_equation_solver: integer to pick the solver used to evaluate sensitivity values in submethods 3 and 0
%       0: Euler method on the sensitivity definition
%       1: Euler method on Jacobian + ODE solver on sensitivity ODE
%       2: ODE solver on sensitivity ODE using explicit Jacobian function
%           Requires Input_files/UP_Jacobian_Function.m
%       NaN: let this function select the best available solver based on 
%           provided files in the './Input_files' folder
%           (Default if no file is provided: 0) 
%       See Sensitivity_solver.m solver_choice for details on each solver.

% List of outputs
%   [S_x_low,S_x_up]: bounds of the sensitivity of successors (at time t_final) to variations of initial states (at time t_init)
%   [S_p_low,S_p_up]: bounds of the sensitivity of successors (at time t_final) to variations of inputs

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function [S_x_low,S_x_up,S_p_low,S_p_up] = Sensitivity_bounds_choice(t_init, t_final, x_low, x_up, p_low, p_up, sensitivity_bounds_method, sensitivity_equation_solver)

% State and input dimensions
n_x = length(x_low);
n_p = length(p_low);

log = logger.get_logger();

%% Hierarchy of the methods to obtain sensitivity bounds

% If sensitivity_bounds_method is undefined (NaN), we follow the order:
% (4-1): Read sensitivity bounds from user-provided files 
%   Requires either Input_files/UP_Sensitivity_Bounds.m 
%   or Input_files/UP_Sensitivity_Signs.m
% (4-2): Interval Arithmetics 
%   Requires Input_files/UP_Jacobian_Function.m
% (4-3): 2nd order Sensitivity 
%   Requires Input_files/UP_Jacobian_Bounds.m 
%   and Input_files/UP_Jacobian_2ndOrder_Bounds.m
% (4-0-2): Sampling and falsification 
%   using sensitivity_equation_solver=2 with explicit Jacobian function
%   Requires Input_files/UP_Jacobian_Function.m
% (4-0-0): Sampling and falsification 
%   using sensitivity_equation_solver=0 with Euler method on the sensitivity

% If sensitivity_bounds_method is provided by the user, we check if the
% requirements associated with the chosen method are fullfiled, and return
% an error message otherwise.
    
METHOD_USER_PROVIDED = 1;
METHOD_INTERVAL_ARITHMETICS = 2;
METHOD_SECOND_ORDER = 3;
METHOD_SAMPLING_FALSIFICATION = 0;

SOLVER_EULER_SENSITIVITY = 0;
SOLVER_EULER_JACOBIAN_ODE_SENSITIVITY = 1;
SOLVER_ODE_SENSITIVITY_EXPLICIT_JACOBIAN = 2;

%% Read sensitivity bounds from a user-provided file
if isnan(sensitivity_bounds_method) || sensitivity_bounds_method == METHOD_USER_PROVIDED
    % Check Input_files/UP_Sensitivity_Bounds.m
    [S_x_low,S_x_up,S_p_low,S_p_up] = UP_Sensitivity_Bounds(t_init,t_final,x_low,x_up,p_low,p_up);
    if ~any(isnan(S_x_low(:))) && ~any(isnan(S_x_up(:))) && ~any(isnan(S_p_low(:))) && ~any(isnan(S_p_up(:)))
        % If all 4 bounds are defined, use them directly
        log.info('(4-1) Extract sensitivity bounds from Input_files/UP_Sensitivity_Bounds.m\n')
        return
    end

    % Check Input_files/UP_Sensitivity_Signs.m
    [S_x_signs,S_p_signs] = UP_Sensitivity_Signs(t_init,t_final,x_low,x_up,p_low,p_up);
    if ~any(isnan(S_x_signs(:))) && ~any(isnan(S_p_signs(:)))
        % If sign matrices are defined, convert them into simple bounds
        % [1,1] if positive, [-1,-1] if negative, [0,0] otherwise
        S_x_low = (S_x_signs>0) - (S_x_signs<0);
        S_x_up = S_x_low;
        S_p_low = (S_p_signs>0) - (S_p_signs<0);
        S_p_up = S_p_low;
        log.info('(4-1) Extract sensitivity signs from Input_files/UP_Sensitivity_Signs.m\n')
        return
    end

    % If both are undefined and user asked for this method
    if sensitivity_bounds_method == METHOD_USER_PROVIDED
        % Return error message to user
        log.error('%s (sensitivity_bounds_method=%d), but:\n%s\n%s', ...
              'Asked to read sensitivity bounds from user-provided files', sensitivity_bounds_method, ...
              ' - sensitivity bounds are not defined in Input_files/UP_Sensitivity_Bounds.m', ...
              ' - sensitivity signs are not defined in Input_files/UP_Sensitivity_Signs.m')
    end
end

%% Interval arithmetics
if isnan(sensitivity_bounds_method) || sensitivity_bounds_method == METHOD_INTERVAL_ARITHMETICS
    % Check Input_files/UP_Jacobian_Bounds.m
    [J_x_low,J_x_up,J_p_low,J_p_up] = UP_Jacobian_Bounds(t_init,t_final,x_low,x_up,p_low,p_up);
    if ~any(isnan(J_x_low(:))) && ~any(isnan(J_x_up(:))) && ~any(isnan(J_p_low(:))) && ~any(isnan(J_p_up(:)))
        if ~any(isinf(J_x_low(:))) && ~any(isinf(J_x_up(:))) && ~any(isinf(J_p_low(:))) && ~any(isinf(J_p_up(:)))
            % The Jacobian bounds are defined and not infinite
           log.info('(4-2) Extract sensitivity bounds with interval arithmetics\n');
           
            % If the user did not explicitely ask for this method 
            % (picked by default due to UP_Jacobian_Bounds.m being defined) 
            if isnan(sensitivity_bounds_method)
                % Inform the user of possible alternative 
                %   if the obtained results are too conservative 
                log.warning('If the final reachability results seem too conservative, consider the tighter but more expensive sensitivity_bounds_method = 3 or 0\n');
            end
            
            % Then over-approximate the sensitivity bounds using interval arithmetics
            [S_x_low,S_x_up,S_p_low,S_p_up] = Sensitivity_interval_arithmetics(t_init,t_final,x_low,x_up,p_low,p_up,J_x_low,J_x_up,J_p_low,J_p_up);
            return
        else
            % If the infinity norm of the Jacobian bounds is inf, then this
            % method cannot be used since the interval matrix exponential
            % will fail.
            log.error('%s (sensitivity_bounds_method=%d), but:\n%s', ...
              'Asked to use interval arithmetics to compute sensitivity bounds', sensitivity_bounds_method, ...
              ' - Jacobian bounds contain inf elements in Input_files/UP_Jacobian_Bounds.m');
        end
    elseif sensitivity_bounds_method == METHOD_INTERVAL_ARITHMETICS
        % If Jacobian bounds are undefined and user asked for this method,
        %   return error message to user
        log.error('%s (sensitivity_bounds_method=%d), but:\n%s', ...
              'Asked to use interval arithmetics to compute sensitivity bounds', sensitivity_bounds_method, ...
              ' - Jacobian bounds are not defined in Input_files/UP_Jacobian_Bounds.m');
    end
end

%% Pick Sensitivity Solver to be used in the last two methods
% 3 solvers of the sensitivity equations are possible
%   0: Euler method on the sensitivity definition
%   1: Euler method on Jacobian + ODE solver on sensitivity ODE
%   2: ODE solver on sensitivity ODE using explicit Jacobian function
%       (requires Input_files/UP_Jacobian_Function.m)
%   See Sensitivity_solver.m solver_choice for details.

% If the user asked to use the sensitivity solver 2 
%   or did not specify which solver to use 
if isnan(sensitivity_equation_solver) || sensitivity_equation_solver == SOLVER_ODE_SENSITIVITY_EXPLICIT_JACOBIAN
    % Check Input_files/UP_Jacobian_Function.m
    [J_x,J_p] = UP_Jacobian_Function(t_init,t_init,x_low,p_low);
    if ~any(isnan(J_x(:))) && ~any(isnan(J_p(:)))
        % If the Jacobian function is defined, use solver 2
        sensitivity_equation_solver = SOLVER_ODE_SENSITIVITY_EXPLICIT_JACOBIAN;
    elseif sensitivity_equation_solver == SOLVER_ODE_SENSITIVITY_EXPLICIT_JACOBIAN
        % If the Jacobian function is not defined 
        %   and user asked for this solver, return an error
        log.error('%s %d but:\n%s', ...
              'Asked to use sensitivity solver', sensitivity_equation_solver, ...
              ' - Jacobian functions are not defined in Input_files/UP_Jacobian_Function.m')
    else
        % Otherwise, pick sensitivity solver 0 by default
        sensitivity_equation_solver = SOLVER_EULER_SENSITIVITY;
    end
end

switch sensitivity_equation_solver
    case SOLVER_EULER_JACOBIAN_ODE_SENSITIVITY
        log.info('(4-x-1) Sensitivity solver used: solve sensitivity ODE with Euler approximation of Jacobians.\n')
    case SOLVER_ODE_SENSITIVITY_EXPLICIT_JACOBIAN
        log.info('(4-x-2) Sensitivity solver used: solve sensitivity ODE with user-provided Jacobian function.\n')
    otherwise
        log.info('(4-x-0) Sensitivity solver used: Euler approximation of sensitivity definition.\n')
end

%% Second-order sensitivity submethod
if isnan(sensitivity_bounds_method) || sensitivity_bounds_method == METHOD_SECOND_ORDER   
    run('Solver_parameters.m'); 

    % Check first-order Jacobian in Input_files/UP_Jacobian_Bounds.m
    [J_x_low, J_x_up, J_p_low, J_p_up] = UP_Jacobian_Bounds(t_init, t_final, x_low, x_up, p_low, p_up); 
    % Check second-order Jacobian in Input_files/UP_Jacobian_2ndOrder_Bounds.m
    [J_xx_low, J_xx_up, J_xp_low, J_xp_up, J_pp_low, J_pp_up, J_px_low, J_px_up] = UP_Jacobian_2ndOrder_Bounds(t_init, t_final, x_low, x_up, p_low, p_up);

    if ~any_is_nan({parameters.sensitivity_samples_per_state_dimension,...
            parameters.sensitivity_samples_per_input_dimension,...
            sensitivity_equation_solver,...
            J_x_low, J_x_up, J_p_low, J_p_up,...
            J_xx_low, J_xx_up, J_xp_low, J_xp_up,...
            J_pp_low, J_pp_up, J_px_low, J_px_up})
        
        log.info('(4-3) Extract sensitivity bounds with second-order sensitivity method.\n');
        
        [S_x_low, S_x_up, S_p_low, S_p_up] = Sensitivity_from_second_order(...
            t_init, t_final,...
            x_low, x_up, p_low, p_up,...
            parameters.sensitivity_samples_per_state_dimension,...
            parameters.sensitivity_samples_per_input_dimension,...
            parameters.sensitivity_equation_solver,...
            J_x_low, J_x_up, J_p_low, J_p_up,...
            J_xx_low, J_xx_up, J_xp_low, J_xp_up,...
            J_pp_low, J_pp_up, J_px_low, J_px_up);
        return;
    elseif sensitivity_bounds_method == METHOD_SECOND_ORDER
        % If some of the required parameters are not defined 
        %   and user asked for this solver, return an error
        log.error('%s (sensitivity_bounds_method=%d), but\n%s\n%s\n%s\n%s\n%s', ...
                'Asked to use second-order sensitivites to compute sensitivity bounds', sensitivity_bounds_method, ...
                'some of the required parameters are undefined amongst:', ...
                ' - sensitivity_samples_per_state_dimension',...
                ' - sensitivity_samples_per_input_dimension',...
                ' - first-order Jacobian in Input_files/UP_Jacobian_Bounds.m',...
                ' - second-order Jacobian in Input_files/UP_Jacobian_2ndOrder_Bounds.m'); 
    end
end

%% Sampling and Falsification (default)
log.info('(4-0) Extract sensitivity bounds with the sampling and falsification method.\n');
log.warning('This submethod does not give a guaranteed OVER-approximation of the reachable set\n');

[S_x_low,S_x_up,S_p_low,S_p_up] = Sensitivity_sampling(t_init, t_final,...
    x_low, x_up, p_low, p_up, sensitivity_equation_solver);

[S_x_low,S_x_up,S_p_low,S_p_up] = Sensitivity_falsification(t_init, t_final,...
    x_low, x_up, p_low, p_up,...
    S_x_low, S_x_up, S_p_low, S_p_up,...
    sensitivity_equation_solver);

end % end Sensitivity_bounds_choice

%% Check if a cell array contains any NaN elements
% List of inputs
%   C: a Cell array
% List of outputs
%   result: true if any element of C is NaN, false otherwise
function result = any_is_nan(C)
    result = false;
    for i = 1:length(C)
       result = result || any(isnan(C{i}), 'all');
    end
end % end any_is_nan
