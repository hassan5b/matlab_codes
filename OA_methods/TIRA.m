%% Hub function serving as interface between the reachability problem definition and the library of over-approximation methods available in folder ./OA_methods
% This function gathers the call of all over-approximation methods, each
% preceded by the tests of their requirements.
% If a specific method is requested by the user (either as the last
% argument of this function call, or defined in Solver_parameters.m), the
% function only checks the requirements for this method and then calls the
% over-approximation function (or returns an error if the requirements are
% not satisfied).
% Otherwise, the function applies the first over-approximation method whose
% requirements are met (in the order detailed below).

% Source:
% P.-J. Meyer, A. Devonport and M. Arcak, "Interval Reachability Analysis: 
% Bounding Trajectories of Uncertain Systems with Boxes for Control and 
% Verification". Springer, 2021. DOI: 10.1007/978-3-030-65110-7

% List of inputs
%   time_vector: time information for the reachability analysis
%       1D vector for discrete-time system: time_vector=t_init
%       2D vector for continuous-time system: time_vector=[t_init,t_final]
%   [x_low,x_up]: interval of initial states (at time t_init)
%   [p_low,p_up]: interval of allowed input values

% Optional input
%   OA_method: integer to request a specific over-approximation method
%       if not provided, checks if it is defined in Solver_parameters.m
%       if not defined in Solver_parameters.m either, the file picks the
%       first method whose requirements are met in the following order:
%       1 - Continuous-time monotonicity (skipped if OA_method=NaN)
%       2 - Continuous-time contraction and growth bound
%       3 - Continuous-time mixed-monotonicity
%       4 - Continuous-time sampled-data mixed-monotonicity
%       5 - Discrete-time mixed-monotonicity
%       6 - Continuous-time and discrete-time Quasi-Monte Carlo
%       7 - Continuous-time and discrete-time Monte Carlo

% List of outputs
%   [succ_low,succ_up]: interval over-approximation of the reachable set
%       after one step for discrete-time system
%       or at time t_final for continuous-time system

% Authors:
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function [succ_low,succ_up] = TIRA(time_vector,x_low,x_up,p_low,p_up,OA_method)
%% Initialization

% State and input dimensions
n_x = length(x_low);
n_p = length(p_low);

% Extract method choices and parameters in the 'parameters' structure
run('Solver_parameters.m');

% Use the over-approximation method choice from Solver_parameters.m
%   if it was not provided as input argument of this function
if nargin < 6
    OA_method = parameters.OA_method;
end

% Create logger
[log, delete_logger] = logger.get_logger();
log.info('Start of TIRA call\n');

%% Extract time information
if length(time_vector) == 1
    % 1D vector: discrete-time system (time_vector = t_init)
    t_init = time_vector;
    t_final = NaN;
elseif length(time_vector) == 2
    % 2D vector: continuous-time system (time_vector = [t_init,t_final])
    t_init = time_vector(1);
    t_final = time_vector(2);
    if t_final == t_init
        % Test if the time range is a single time instant
        log.info('The initial and final times are equals: the reachable set is the interval of initial states\n')
        succ_low = x_low;
        succ_up = x_up;
        delete_logger();
        return
    elseif t_final < t_init
        % Test if the time range is not properly defined
        log.error('Inconsistent time range provided (t_final < t_init)')
    end
else
    log.error('time_vector can only contain 1 (for discrete-time system) or 2 elements (for continuous-time system)')
end

%% Check if the system is defined
system_eval = System_description(t_init,x_low,p_low);
assert(~any(isnan(system_eval)),'The considered system needs to be defined in System_description.m')

%% Check consistency of user-provided intervals
assert(all(x_low<=x_up),'Inconsistent interval of initial states provided: need x_low <= x_up')
assert(all(p_low<=p_up),'Inconsistent interval of inputs provided: need p_low <= p_up')

%% Check if over-approximation of the reachable set can be skipped
if all(isequal(x_low,x_up)) && all(isequal(p_low,p_up))
    % Both interval contain a single point: simply compute the successor
    log.info('Provided intervals are singletons: exact reachable set is a singleton too\n')
    if isnan(t_final)
        % Discrete-time successor
        succ_low = System_description(t_init,x_low,p_low);
        succ_up = succ_low;
    else
        % Continuous-time successor
        ode_solver = parameters.ode_solver;
        ode_options = parameters.ode_options;
        [~,x_traj] = ode_solver(@(t,x) System_description(t,x,p_low), [t_init, t_final], x_low, ode_options);
        succ_low = x_traj(end,:)';
        succ_up = succ_low;
    end
    delete_logger();
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous-time methods %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - Continuous-time monotonicity
% This choice will check if the system is monotone, then apply an
% over-approximation method specific to monotone systems.
% This method is only checked if user specifically asks for it in OA_method
% (by construction, it is encompassed in the less scalable OA_method=3)

% *Requirements:
% - Continuous-time system
% - One of the following 2 user-provided items:
%     Signs of the Jacobian matrices in ./Input_files/UP_Jacobian_Signs.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m

if ~isnan(t_final) && OA_method == 1
    tic_CTM = tic;
    log.info('(1) Checking conditions for continuous-time monotonicity method...\n')
    
    % Check Input_files/UP_Jacobian_Signs.m
    [J_x_signs,J_p_signs] = UP_Jacobian_Signs(t_init,t_final,x_low,x_up,p_low,p_up);
    
    if any(isnan(J_x_signs(~logical(eye(n_x))))) || any(isnan(J_p_signs(:)))
        % If the Jacobian signs are not defined,
        % Check Input_files/UP_Jacobian_Bounds.m
        [J_x_low,J_x_up,J_p_low,J_p_up] = UP_Jacobian_Bounds(t_init,t_final,x_low,x_up,p_low,p_up);
        if any(isnan(J_x_low(~logical(eye(n_x))))) || any(isnan(J_x_up(~logical(eye(n_x))))) || any(isnan(J_p_low(:))) || any(isnan(J_p_up(:)))
            % If the Jacobian bounds are not defined
            % Return an error since this method was specifically requested
            log.error('%s (OA_method=%d), but:\n%s\n%s', ...
                'Continuous-time monotonicitiy approach selected', OA_method, ...
                ' - Jacobian signs are not defined in Input_files/UP_Jacobian_Signs.m', ...
                ' - Jacobian bounds are not defined in Input_files/UP_Jacobian_Bounds.m')
        else
            log.info('Read Jacobian bounds from Input_files/UP_Jacobian_Bounds.m\n')
            
            % Check consistency of user-provided intervals
            assert(all(J_x_low(~logical(eye(n_x)))<=J_x_up(~logical(eye(n_x)))),'Inconsistent bounds of the state Jacobian provided: need J_x_low <= J_x_up')
            assert(all(J_p_low(:)<=J_p_up(:)),'Inconsistent bounds of the input Jacobian provided: need J_p_low <= J_p_up')
            
            % Check sign-stability of the Jacobians bounds
            if all(J_x_low(~logical(eye(n_x)))>=0 | J_x_up(~logical(eye(n_x)))<=0) && all(J_p_low(:)>=0 | J_p_up(:)<=0)
                % If sign-stable, extract the corresponding signs
                J_x_signs = (J_x_up>0) - (J_x_low<0);
                J_p_signs = (J_p_up>0) - (J_p_low<0);
                log.info('Jacobian matrices are sign-stable.\n')
            else
                % If not sign-stable, return an error 
                %   since this method was specifically requested
                log.error('%s (OA_method=%d), but:\n%s', ...
                    'Continuous-time monotonicitiy approach selected', OA_method, ...
                    ' - Jacobian bounds defined in Input_files/UP_Jacobian_Bounds.m are not sign-stable')
            end
        end
    else
        log.info('Read Jacobian signs from Input_files/UP_Jacobian_Signs.m\n')
    end
    
    
    % Check the monotonicity condition on the Jacobian signs
    
    % Continuous-time monotonicity has no condition on diag(J_x_signs)
    % To avoid leaving NaN diagonal elements, replace diagonal by zeros
    J_x_signs_no_nan = J_x_signs;
    J_x_signs_no_nan(1:n_x+1:end) = 0;

    % Check if Jacobian sign matrices correspond to a monotone system
    [is_monotone,partial_order_x,partial_order_p] = Monotonicity_sign_conditions(J_x_signs_no_nan,J_p_signs);
    if is_monotone
        log.info('Jacobian signs satisfy the monotonicity condition.\n')
        log.info('Call of the continuous-time monotonicity over-approximation method ...\n')
        % Call over-approximation method
        [succ_low,succ_up] = OA_1_CT_Monotonicity(t_init,t_final,x_low,x_up,p_low,p_up,partial_order_x,partial_order_p);
        toc_CTM = toc(tic_CTM);
        log.runtime('Checking monotonicity and calling the continuous-time monotonicity method took %f seconds\n', toc_CTM);
        delete_logger();
        return
    else
        % Return an error since this method was specifically requested
        log.error('%s (OA_method=%d), but:\n%s', ...
            'Continuous-time monotonicitiy approach selected', OA_method, ...
            ' - provided Jacobian signs or bounds do not correspond to a monotone system')
    end
end

%% 2 - Continuous-time contraction and growth bound
% This method relies on using or defining a function that bounds the
% growth or contraction on each dimension of the system.

% *Requirements:
% - Continuous-time system
% - One of the following 3 user-provided items:
%     Growth bound function in ./Input_files/UP_Growth_Bound_Function.m
%     Contraction matrix or scalar in ./Input_files/UP_Contraction_Matrix.m
%     Bounds of the state Jacobian in ./Input_files/UP_Jacobian_Bounds.m
% *The last 2 items also require dynamics with additive input:
%   System_description(t,x,p)=System_description(t,x,0)+p
% which is automatically checked in ./Utilities/Growth_bound_choice.m

if ~isnan(t_final) && (isnan(OA_method) || OA_method == 2)
    tic_GB = tic;
    log.info('(2) Checking conditions for continuous-time method based on contraction or growth bounds...\n')
    
    % Check the requirements of the growth bound method
    %   and if satisfied, extract a growth bound function handle
    [GB_handle,bool_skip_CTGB] = Growth_bound_choice(t_init,t_final,x_low,x_up,p_low,p_up);
    
    if ~bool_skip_CTGB
        % If a growth bound handle was successfully obtained
        log.info('Call of the contraction/growth bound over-approximation method ...\n')
        % Call over-approximation method
        [succ_low,succ_up] = OA_2_CT_Contraction_growth_Bound(t_init,t_final,x_low,x_up,p_low,p_up,GB_handle);
        toc_GB = toc(tic_GB);
        log.runtime('Checking growth-bound conditions and calling the method took %f seconds\n', toc_GB);
        delete_logger();
        return
    elseif OA_method == 2
        % Otherwise return an error if this specific method was requested
        if bool_skip_CTGB == 1
            log.error('%s (OA_method=%d), but:\n%s\n%s', ...
                'Continuous-time contraction/growth bound approach selected', OA_method, ...
                ' - growth bound function is not defined in Input_files/UP_Growth_Bound_Function.m', ...
                ' - the dynamics do not have additive input: dx/dt=f(t,x,p)~=f(t,x,0)+p')
        else
            log.error('%s (OA_method=%d), but:\n%s\n%s\n%s', ...
                'Continuous-time contraction/growth bound approach selected', OA_method, ...
                ' - growth bound function is not defined in Input_files/UP_Growth_Bound_Function.m', ...
                ' - contraction matrix (or scalar) is not defined in Input_files/UP_Contraction_Matrix.m', ...
                ' - Jacobian bounds are not defined in Input_files/UP_Jacobian_Bounds.m')
        end
    else
        % Otherwise, move on to checking the next method
        toc_GB = toc(tic_GB);
        log.runtime('Checking growth-bound requirements took %f seconds\n', toc_GB);
    end
end

%% 3 - Continuous-time mixed-monotonicity
% This method is an extension of method 1 to any system with bounded
% Jacobian matrices.

% *Requirements:
% - Continuous-time system
% - One of the following 2 user-provided items:
%     Signs of the Jacobian matrices in ./Input_files/UP_Jacobian_Signs.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m

bool_skip_CTMM = 0;
if ~isnan(t_final) && (isnan(OA_method) || OA_method == 3)
    tic_CTMM = tic;
    log.info('(3) Checking conditions for continuous-time mixed-monotonicity method...\n')
    
    % Check Input_files/UP_Jacobian_Bounds.m
    [J_x_low,J_x_up,J_p_low,J_p_up] = UP_Jacobian_Bounds(t_init,t_final,x_low,x_up,p_low,p_up);
    if any(isnan(J_x_low(~logical(eye(n_x))))) || any(isnan(J_x_up(~logical(eye(n_x))))) || any(isnan(J_p_low(:))) || any(isnan(J_p_up(:)))
        % If the Jacobian bounds are not defined,
        % Check Input_files/UP_Jacobian_Signs.m
        [J_x_signs,J_p_signs] = UP_Jacobian_Signs(t_init,t_final,x_low,x_up,p_low,p_up);
        if any(isnan(J_x_signs(~logical(eye(n_x))))) || any(isnan(J_p_signs(:)))
            % If the Jacobian signs are not defined either
            if OA_method == 3
                % Return an error if this method was specifically requested
                log.error('%s (OA_method=%d), but:\n%s\n%s', ...
                    'Continuous-time mixed-monotonicitiy approach selected', OA_method, ...
                    ' - Jacobian bounds are not defined in Input_files/UP_Jacobian_Bounds.m', ...
                    ' - Jacobian signs are not defined in Input_files/UP_Jacobian_Signs.m')
            else
                log.info('Neither Jacobian bounds nor signs are provided.\n')
                % Flag to go to the next method if OA_method was undefined
                bool_skip_CTMM = 1;
            end
        else
            % If sign matrices are defined, convert them into simple bounds
            % [1,1] if positive, [-1,-1] if negative, [0,0] otherwise
            J_x_low = (J_x_signs>0) - (J_x_signs<0);
            J_x_up = J_x_low;
            J_p_low = (J_p_signs>0) - (J_p_signs<0);
            J_p_up = J_p_low;
            log.info('Read Jacobian signs from Input_files/UP_Jacobian_Signs.m\n')
        end
    else
        log.info('Read Jacobian bounds from Input_files/UP_Jacobian_Bounds.m\n')
    end
    
    % Only apply the over-approximation method if we have Jacobian bounds
    if ~bool_skip_CTMM
        log.info('Call of the continuous-time mixed-monotonicity over-approximation method ...\n')
        % Call over-approximation method
        [succ_low,succ_up] = OA_3_CT_Mixed_Monotonicity(t_init,t_final,x_low,x_up,p_low,p_up,J_x_low,J_x_up,J_p_low,J_p_up);
        toc_CTMM = toc(tic_CTMM);
        log.runtime('Checking conditions and calling continuous-time mixed-monotonicity method took %f seconds\n', toc_CTMM);
        delete_logger();
        return
    else
        % Otherwise, move on to checking the next method
        toc_CTMM = toc(tic_CTMM);
        log.runtime('Checking conditions for continuous-time mixed-monotonicity method took %f seconds\n', toc_CTMM);
    end
end

%% 4 - Continuous-time sampled-data mixed-monotonicity
% This method applies the discrete-time mixed-monotonicity method (5)
% to the sampled-data version of a continuous-time system.
% The requirements on the Jacobian of the discrete-time system thus
% translate into requirements on the sensitivity matrices of the 
% continuous-time system.

% *Requirements:
% - Continuous-time system
% - Input p constant over the considered time range [t_init,t_final].
% - Optional user-provided inputs amongst the following items:
%     Signs of the sensitivity matrices in ./Input_files/UP_Sensitivity_Signs.m
%     Bounds of the sensitivity matrices in ./Input_files/UP_Sensitivity_Bounds.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m
%     Bounds of the 2nd order Jacobian matrices in ./Input_files/UP_Jacobian_2ndOrder_Bounds.m
%     Jacobian functions in ./Input_files/UP_Jacobian_Function.m
% *If none of these optional input is provided, the method can still
% be called using default submethod 4-0 "Sampling and falsification".

if ~isnan(t_final) && (isnan(OA_method) || OA_method == 4)
    tic_SDMM = tic;
    log.info('(4) Checking conditions for sampled-data mixed-monotonicity method based...\n')
    log.warning('This method over-approximates the reachable set only for constant input functions.\n')
    
    % Extract or compute sensitivity bounds
    % (submethod used depends on the available user-provided inputs)
    [S_x_low,S_x_up,S_p_low,S_p_up] = Sensitivity_bounds_choice(...
        t_init, t_final, x_low, x_up, p_low, p_up,...
        parameters.sensitivity_bounds_method,...
        parameters.sensitivity_equation_solver);
    
    % Call over-approximation method
    log.info('Call of the sampled-data mixed-monotonicity over-approximation method ...\n')
    [succ_low,succ_up] = OA_4_CT_Sampled_data_MM(t_init,t_final,x_low,x_up,p_low,p_up,S_x_low,S_x_up,S_p_low,S_p_up);
    
    toc_SDMM = toc(tic_SDMM);
    log.runtime('Checking conditions and calling sampled-data mixed-monotonicity method took %f seconds\n', toc_SDMM);
    delete_logger();
    return;
    
    % *Reasons why we do not check other continuous-time methods down this
    % file after this method:
    
    % - If the user did not request specific submethod or sensitivity 
    % solver for this method (in Solver_parameters.m), then this method can
    % always be used with the "Sampling and falsification" submethod (4-0) 
    % with the default Euler sensitivity solver (4-0-0).
    
    % - If the user requested specific submethod or sensitivity solver for 
    % this method in Solver_parameters.m (even if they did not specifically
    % ask for this main method (OA_method=NaN)), then the code above only 
    % check the requested submethod: applies it if its requirements are 
    % met; returns an error if they are not met.
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete-time methods %
%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5 - Discrete-time mixed-monotonicity
% This method is similar to method 3 (i.e. applicable to any system with
% bounded Jacobian matrices) but for discrete-time systems.

% *Requirements:
% - Discrete-time system
% - One of the following 2 user-provided items:
%     Signs of the Jacobian matrices in ./Input_files/UP_Jacobian_Signs.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m

bool_skip_DTMM = 0;
if isnan(t_final) && (isnan(OA_method) || OA_method == 5)
    tic_DTMM = tic;
    log.info('(5) Checking conditions for discrete-time mixed-monotonicity method...\n')
    
    % Check Input_files/UP_Jacobian_Bounds.m
    [J_x_low,J_x_up,J_p_low,J_p_up] = UP_Jacobian_Bounds(t_init,t_final,x_low,x_up,p_low,p_up);
    if any(isnan(J_x_low(:))) || any(isnan(J_x_up(:))) || any(isnan(J_p_low(:))) || any(isnan(J_p_up(:)))
        % If the Jacobian bounds are not defined,
        % Check Input_files/UP_Jacobian_Signs.m
        [J_x_signs,J_p_signs] = UP_Jacobian_Signs(t_init,t_final,x_low,x_up,p_low,p_up);
        if any(isnan(J_x_signs(:))) || any(isnan(J_p_signs(:)))
            % If the Jacobian signs are not defined either
            if OA_method == 5
                % Return an error if this method was specifically requested
                log.error('%s (OA_method=%d), but:\n%s\n%s', ...
                    'Discrete-time mixed-monotonicitiy approach selected', OA_method, ...
                    ' - Jacobian bounds are not defined in Input_files/UP_Jacobian_Bounds.m', ...
                    ' - Jacobian signs are not defined in Input_files/UP_Jacobian_Signs.m')
            else
                log.info('Neither Jacobian bounds nor signs are provided.\n')
                % Flag to go to the next method if OA_method was undefined
                bool_skip_DTMM = 1;
            end
        else
            % If sign matrices are defined, convert them into simple bounds
            % [1,1] if positive, [-1,-1] if negative, [0,0] otherwise
            J_x_low = (J_x_signs>0) - (J_x_signs<0);
            J_x_up = J_x_low;
            J_p_low = (J_p_signs>0) - (J_p_signs<0);
            J_p_up = J_p_low;
            log.info('Read Jacobian signs from Input_files/UP_Jacobian_Signs.m\n')
        end
    else
        log.info('Read Jacobian bounds from Input_files/UP_Jacobian_Bounds.m\n')
    end
    
    % Only apply the over-approximation method if we have Jacobian bounds
    if ~bool_skip_DTMM
        log.info('Call of the discrete-time mixed-monotonicity over-approximation method ...\n')
        % Call over-approximation method
        [succ_low,succ_up] = OA_5_DT_Mixed_Monotonicity(t_init,x_low,x_up,p_low,p_up,J_x_low,J_x_up,J_p_low,J_p_up);
        toc_DTMM = toc(tic_DTMM);
        log.runtime('Checking conditions and calling discrete-time mixed-monotonicity method took %f seconds\n', toc_DTMM);
        delete_logger();
        return
    else
        % Otherwise, move on to checking the next method
        toc_DTMM = toc(tic_DTMM);
        log.runtime('Checking conditions for discrete-time mixed-monotonicity method took %f seconds\n', toc_DTMM);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete and Continuous-time methods %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6 - Quasi-Monte Carlo
% This method simulates the system over a grid of points and combines
% it with an error-bound function to over-approximate the reachable set.

% *Requirements for discrete-time systems:
% - Error-bound function in ./Input_files/UP_Error_Bound_Function.m

% *Requirements for continuous-time systems:
% - One of the following 4 user-provided items:
%     Error-bound function in ./Input_files/UP_Error_Bound_Function.m
%     Lipschitz constant in ./Input_files/UP_Lipschitz_Constant.m
%     Contraction matrix or scalar in ./Input_files/UP_Contraction_Matrix.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m
% *The last 2 items also require dynamics with additive input:
%   System_description(t,x,p)=System_description(t,x,0)+p
% which is automatically checked in ./Utilities/Error_bound_choice.m

if isnan(OA_method) || OA_method == 6
    tic_QMC = tic;
    log.info('(6) Checking conditions for quasi-Monte Carlo method...\n')
    
    % Check the requirements of the error bound method
    %   and if satisfied, extract an error bound function handle
    [error_bound_handle, bool_skip_QMC] = Error_bound_choice(time_vector, x_low, x_up, p_low, p_up, parameters.error_bound_method);
    
    if ~bool_skip_QMC
        % If an error bound handle was successfully obtained
        log.info('Call of the quasi-Monte Carlo over-approximation method ...\n')
        % Call over-approximation method
        [succ_low, succ_up] = OA_6_DT_CT_Quasi_Monte_Carlo(time_vector,...
            x_low, x_up, p_low, p_up,...
            parameters.QMC_samples_per_state_dimension,...
            parameters.QMC_samples_per_input_dimension,...
            error_bound_handle);
        toc_QMC = toc(tic_QMC);
        log.runtime('Checking error bound conditions and calling the quasi-Monte Carlo method took %f seconds\n', toc_QMC);
        delete_logger();
        return
    elseif OA_method == 6
        % Otherwise return an error if this specific method was requested
        if bool_skip_QMC == 1
            log.error('%s (OA_method=%d), but:\n%s',...
                'Quasi-Monte Carlo approach selected', OA_method,...
                ' - Discrete-time system, but no user-provided error bound function');
        elseif bool_skip_QMC == 2
            log.error('%s (OA_method=%d), but:\n%s\n%s\n%s',...
                'Quasi-Monte Carlo approach selected', OA_method,...
                ' - Contraction matrix submethod fails since input additivity requirement (that dx/dt = f(t,x,p) = f(t,x,0) + p) not satisfied',...
                ' - Jacobian bounds submethod fails since input additivity requirement (that dx/dt = f(t,x,p) = f(t,x,0) + p) not satisfied',...
                ' - Lipschitz constant submethod fails since no constant provided in Input_files/UP_Lipschitz_Constant.m');
        elseif bool_skip_QMC == 3
            log.error('%s (OA_method=%d), but:\n%s\n%s\n%s',...
                'Quasi-Monte Carlo approach selected', OA_method,...
                ' - Contraction matrix submethod fails since no matrix provided in Input_files/UP_Contraction_Matrix.m',...
                ' - Jacobian bounds submethod fails since no bound provided in Input_files/UP_Jacobian_Bounds.m',...
                ' - Lipschitz constant submethod fails since no constant provided in Input_files/UP_Lipschitz_Constant.m');
        end
    else
        % Otherwise, move on to checking the next method
        toc_QMC = toc(tic_QMC);
        log.runtime('Checking quasi-Monte Carlo requirements took %f seconds\n', toc_QMC);
    end
end

%% 7 - Monte Carlo
% This method simulates the system over uniformly sampled points to 
% generate an approximation interval with a desired probabilistic 
% accuracy and confidence.
% *The resulting interval approximation is not guaranteed to be 
% an OVER-approximation of the reachable set.

% *This method has no requirement and is always applicable to both
% discrete-time and continuous-time systems.

if isnan(OA_method) || OA_method == 7
    tic_MonteCarlo = tic;
    log.info('(7) Call of the Monte Carlo approximation method ...\n');
    log.warning('This method does not give a guaranteed OVER-approximation of the reachable set\n');
    % Call approximation method.
    [succ_low, succ_up] = OA_7_DT_CT_Monte_Carlo(time_vector,...
        x_low, x_up, p_low, p_up,...
        parameters.epsilon, parameters.delta);
    toc_MonteCarlo = toc(tic_MonteCarlo);
    log.runtime('Call of the Monte Carlo method took %f seconds\n', toc_MonteCarlo);
    delete_logger();
    return
end

%% Error message if we reach the end of this file (this should not be possible)
log.error('None of the over-approximation methods could be applied.')
