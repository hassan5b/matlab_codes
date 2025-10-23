%% Extract or create an error bound function handle
% The user can either ask for a specific method, or let the function pick
% a method based on which user-proviced files are available in the folder
% ./Input_files

% Source paper:
% P.-J. Meyer, A. Devonport and M. Arcak, "Sampling-Based Methods", 
% In: "Interval Reachability Analysis: Bounding Trajectories of Uncertain 
% Systems with Boxes for Control and Verification". Springer, 2021
% DOI: 10.1007/978-3-030-65110-7_7

% List of Inputs:
%   time_vector: time information for the reachability analysis
%       1D vector for discrete-time system: time_vector = t_init
%       2D vector for continuous-time system: time_vector=[t_init, t_final]
%       t_init: initial time
%       t_final: if present, time at which reachable set is approximated.
%           For discrete-time systems, t_final will be one step from t_init
%   [x_low, x_up]: interval of initial states (at time t_init)
%   [p_low, p_up]: interval of allowed input values

% Optional input (read it from Solver_parameters.m if not provided in the call)
%   error_bound_method: integer to pick the method to obtain the error bound
%       1: Read it from user-provided file 
%           Requires Input_files/UP_Error_Bound_Function.m
%       2: Create it from user-provided contraction matrix or scalar
%           Requires Input_files/UP_Contraction_Matrix.m
%       3: Create it from user-provided Jacobian bounds 
%           Requires Input_files/UP_Jacobian_Bounds.m)
%       4: Create it from user-provided Lipschitz constant
%           Requires Input_files/UP_Lipschitz_Constant.m
%       NaN: let this function select the best available method (in the 
%           above order) based on provided files in './Input_files' folder

% List of outputs:
%   error_bound_handle: function handle of the error bound
%       inputs: 
%           time_vector, 
%           1/2 * infinity norm of spacing between state sample points,
%           1/2 * infinity norm of spacing between input sample points
%       outputs: upper bound on error
%   skip_quasi_monte_carlo: flag on why this method cannot be applied
%       0: No errors, some method was run.
%       1: No method was run: discrete-time, no user-provided error bound
%           function.
%       2: No method was run: continuous-time, dynamics are not additive
%           and no user-provided inputs necessary for other methods
%           were provided.
%       3: No method was run: continuous-time, dynamics additive but no 
%           user-provided inputs necessary for any method was provided.

% Authors:
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [error_bound_handle, skip_quasi_monte_carlo] = Error_bound_choice(time_vector, x_low, x_up, p_low, p_up, error_bound_method)

% State and input dimensions
n_x = length(x_low);
n_p = length(p_low);

log = logger.get_logger();

% Extract from Solver_parameters.m the method choice
%   if not provided as input argument of this function
if nargin < 6
    run('Solver_parameters.m');
    error_bound_method = parameters.error_bound_method;
end

METHOD_USER_PROVIDED = 1;
METHOD_CONTRACTION_MATRIX = 2;
METHOD_JACOBIAN_BOUNDS = 3;
METHOD_LIPSCHITZ_CONSTANT = 4;

skip_quasi_monte_carlo = 0;

% Verify time vector size.
assert(length(time_vector) == 1 || length(time_vector) == 2,...
    'time_vector must be of length 1 for discrete-time, and of length 2 for continuous-time');

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

% Check that error_bound_method is one of the available methods.
if ~isnan(error_bound_method) && ~ismember(error_bound_method,...
    [METHOD_USER_PROVIDED, METHOD_CONTRACTION_MATRIX,...
    METHOD_JACOBIAN_BOUNDS, METHOD_LIPSCHITZ_CONSTANT])
    log.error('Specified error_bound_method is invalid');
end

%% Use directly the user-provided error bound function
if isnan(error_bound_method) || error_bound_method == METHOD_USER_PROVIDED
    % Check Input_files/UP_Error_Bound_Function.m
    error_bound_handle = UP_Error_Bound_Function(time_vector, x_low, x_up, p_low, p_up);
    if ~any(isnan(error_bound_handle(time_vector, ones(n_x, 1), ones(n_p, 1))))
        % If the function is defined, use it directly
        log.info('(6-1) Extract error bound function handle from Input_files/UP_Error_Bound_Function.m\n')
        return
    end
    
    % If error_bound_handle undefined and user asked for this method
    if error_bound_method == METHOD_USER_PROVIDED
        % Return error message to user
        log.error('%s%d, but:\n%s', ...
            'Asked for error_bound_method=', error_bound_method,...
            ' - error bound function is not defined in Input_files/UP_Error_Bound_Function.m');
    end
end


%% Check if the system has additive input for 2 of the submethods below
% For the next two submethods (6-2) with contraction matrix and (6-3) with
% Jacobian bounds, the dynamics need to be defined with additive input.
% (use tolerance eps(2)=2^(-51) to account for numerical errors)

additive_input = true;
% This test is only relevant to submethods for continuous-time systems
if is_continuous(time_vector)
    t_final = time_vector(2);
    % Check if the dynamics are NOT with additive input
    if n_p ~= n_x || ...
       any(abs((System_description(t_final,x_low,p_up)-System_description(t_final,x_low,zeros(n_p,1))) - p_up) > eps(2) | ...
           abs((System_description(t_final,x_low,p_low)-System_description(t_final,x_low,zeros(n_p,1))) - p_low) > eps(2) | ...
           abs((System_description(t_final,x_up,p_up)-System_description(t_final,x_up,zeros(n_p,1))) - p_up) > eps(2) | ...
           abs((System_description(t_final,x_up,p_low)-System_description(t_final,x_up,zeros(n_p,1))) - p_low) > eps(2))
       % Raise a flag to be checked later
       additive_input = false;
    end
end

%% Create error-bound function from user-provided contraction matrix
if isnan(error_bound_method) || error_bound_method == METHOD_CONTRACTION_MATRIX
    % This method is only relevant to continuous-time systems
    if is_continuous(time_vector)
        % And it requires the system dynamics to be with additive input
        if additive_input
            t_init = time_vector(1);
            t_final = time_vector(2);

            % Check Input_files/UP_Contraction_Matrix.m
            C = UP_Contraction_Matrix(t_init, t_final, x_low, x_up, p_low, p_up);
            if ~any(isnan(C(:))) && (isequal(size(C),[n_x n_x]) || isscalar(C))
                log.info('(6-2) Extract error bound function from contraction matrix in Input_files/UP_Contraction_Matrix.m\n')
                
                % Check consistency of user-provided matrix
                assert(all(C(~logical(eye(n_x)))>=0),'The contraction matrix needs non-negative off-diagonal elements.')
                
                % Create error-bound function from the contraction matrix
                error_bound_handle = error_bound_from_contraction_matrix(C);
                return;
            elseif error_bound_method == METHOD_CONTRACTION_MATRIX
                % If the contraction matrix is undefined and the user asked
                % for this specific submethod, return an error message
                log.error('%s%d, but:\n%s',...
                  'Asked for error_bound_method=', error_bound_method,...
                  ' - contraction matrix (or scalar) is not defined in Input_files/UP_Contraction_Matrix.m');
            end
        elseif error_bound_method == METHOD_CONTRACTION_MATRIX
            % If the system does not have additive input and the user asked
            % for this specific submethod, return an error message
            log.error('%s%d, but:\n%s',...
              'Asked for error_bound_method=', error_bound_method,...
              ' - system dynamics do not have additive input');
        end
    elseif error_bound_method == METHOD_CONTRACTION_MATRIX
        % If the system is not continuous-time and the user asked for this
        % specific submethod, return an error message
        log.error('%s%d, but:\n%s',...
            'Asked for error_bound_method=', error_bound_method,...
            ' - system is not continuous-time');
    end
end

%% Create error-bound function from user-provided Jacobians bounds
% To be usable in the over-approximation using the error bound associated 
% with this submethod, the bounds on the state Jacobian below need to be
% valid not only for initial states in [x_low,x_up], but also for any state
% reachable over time range [t_init,t_final] by the system with constant
% input equal to the center of the input interval [p_low,p_up]

if isnan(error_bound_method) || error_bound_method == METHOD_JACOBIAN_BOUNDS
    % This method is only relevant to continuous-time systems
    if is_continuous(time_vector)
        % And it requires the system dynamics to be with additive input
        if additive_input
            t_init = time_vector(1);
            t_final = time_vector(2);

            % Check Input_files/UP_Jacobian_Bounds.m
            input_center = (p_low + p_up) / 2;
            [J_x_low, J_x_up, ~, ~] = UP_Jacobian_Bounds(t_init, t_final, x_low, x_up, input_center, input_center);
            if ~any(isnan(J_x_low(:))) && ~any(isnan(J_x_up(:)))
                log.info('(6-3) Extract error bound function from state Jacobian bounds in Input_files/UP_Jacobian_Bounds.m\n')

                % Check consistency of user-provided interval
                assert(all(J_x_low(:) <= J_x_up(:)),'Inconsistent bounds of the state Jacobian provided: need J_x_low <= J_x_up')

                % If state Jacobian bounds are defined, 
                % convert them into a contraction matrix.
                C = max(abs(J_x_low),abs(J_x_up));   % Off-diagonal elements: upper bound on the absolute values
                C(1:n_x+1:end) = diag(J_x_up);       % Diagonal elements: take Jacobian upper bound directly

                % Use this matrix to create an error bound function.  
                error_bound_handle = error_bound_from_contraction_matrix(C);
                return
            elseif error_bound_method == METHOD_JACOBIAN_BOUNDS
                % If the Jacobian bounds are undefined and the user asked
                % for this specific submethod, return an error message
                log.error('%s%d, but:\n%s',...
                  'Asked for error_bound_method=', error_bound_method,...
                  ' - Jacobian bounds are not defined in Input_files/UP_Jacobian_Bounds.m');
            end
        elseif error_bound_method == METHOD_JACOBIAN_BOUNDS
            % If the system does not have additive input and the user asked
            % for this specific submethod, return an error message
            log.error('%s%d, but:\n%s',...
              'Asked for error_bound_method=', error_bound_method,...
              ' - system dynamics do not have additive input');
        end
    elseif error_bound_method == METHOD_JACOBIAN_BOUNDS
        % If the system is not continuous-time and the user asked for this
        % specific submethod, return an error message
        log.error('%s%d, but:\n%s',...
            'Asked for error_bound_method=', error_bound_method,...
            ' - system is not continuous-time');
    end
end

%% Create error-bound function from user-provided Lipschitz constant
if isnan(error_bound_method) || error_bound_method == METHOD_LIPSCHITZ_CONSTANT
    % This method is only relevant to continuous-time systems
    if is_continuous(time_vector)
        % Check Input_files/UP_Lipschitz_Constant.m
        L = UP_Lipschitz_Constant(time_vector(1), time_vector(2), x_low, x_up, p_low, p_up);
        if ~isnan(L)
            log.info('(6-4) Extract error bound function from Lipschitz constant in Input_files/UP_Lipschitz_Constant.m\n')

            % Check consistency of user-provided constant
            assert(isscalar(L), 'Lipschitz constant must be a scalar');
            assert(L >= 0, 'Lipschitz constant must be nonnegative');

            % Use this constant to create an error bound function.  
            error_bound_handle = @(time_vector, state_error, input_error)...
                exp(L*(time_vector(2) - time_vector(1))) * ...
                (state_error + input_error);
            return
        elseif error_bound_method == METHOD_LIPSCHITZ_CONSTANT
            % If the Lipschitz constant is undefined and the user asked
            % for this specific submethod, return an error message
            log.error('%s%d, but:n%s',...
                'Asked for error_bound_method=', error_bound_method,...
                ' - Lipschitz constant is not defined in Input_files/UP_Lipschitz_Constant.m');
        end
    elseif error_bound_method == METHOD_LIPSCHITZ_CONSTANT
        % If the system is not continuous-time and the user asked for this
        % specific submethod, return an error message
        log.error('%s%d, but:\n%s',...
            'Asked for error_bound_method=', error_bound_method,...
            ' - system is not continuous-time');
    end
end

%% If error_bound_method was not provided and none of the above method could be applied
% raise a flag that the error bound method cannot be used

% Define the interger flag on why this method cannot be applied
if is_discrete(time_vector)
    % 1: No method was run: discrete-time, no user-provided error bound 
    %   function.
    skip_quasi_monte_carlo = 1;
elseif ~additive_input
    % 2: No method was run: continuous-time, dynamics are not additive and 
    %   no user-provided inputs necessary for other methods were provided.
    skip_quasi_monte_carlo = 2;
else
    % 3: No method was run: continuous-time, dynamics additive but no  
    %   user-provided inputs necessary for any method was provided.
    skip_quasi_monte_carlo = 3;
end

log.info('The requirements for the quasi-Monte Carlo method are not satisfied\n');
error_bound_handle = NaN;
end % end Error_bound_choice


%% Obtaining an error-bound function from a contraction matrix
function error_bound_handle = error_bound_from_contraction_matrix(C)     
    % If not a scalar, use the infinity-matrix measure to convert C to a scalar.
    if ~isscalar(C)
        assert(size(C, 1) == size(C, 2),...
            'Contraction matrix must be size of state vector by size of state vector');
        C = infinity_matrix_measure(C);
    end
    error_bound_handle = error_bound_from_contraction_scalar(C);
end % end error_bound_from_contraction_matrix

%% Obtaining an error-bound function from a contraction scalar
function error_bound_handle = error_bound_from_contraction_scalar(c)
    error_bound_handle = @(time_vector, state_error, input_error)...
        exp(c * (time_vector(2) - time_vector(1))) * state_error...
        + ((exp(c * (time_vector(2) - time_vector(1))) - 1)/c) * input_error;
end % end error_bound_from_contraction_scalar

%% Infinity-matrix measure of a matrix
function measure = infinity_matrix_measure(A)
    measure = max(sum(abs(A) - abs(diag(diag(A))) + diag(diag(A)), 2));
end % end infinity_matrix_measure

%% Test if the time vector corresponds to a continuous-time system
function out = is_continuous(time_vector)
    if length(time_vector) == 2
        out = true;
    else
        out = false;
    end
end % end is_continuous

%% Test if the time vector corresponds to a discrete-time system
function out = is_discrete(time_vector)
    if length(time_vector) == 1
        out = true;
    else
        out = false;
    end
end % end is_discrete
