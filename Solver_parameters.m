%% Choice of over-approximation methods, submethods and definition of their internal parameters
% All variables defined here are stored in the structure 'parameter', which
% is partially read in the relevant functions.
% The user can modify any of the variables defined in this file:
%  - either to ask for a specific over-approximation method or submethod
%  - or to change the internal parameters of these methods.

% Source:
% P.-J. Meyer, A. Devonport and M. Arcak, "Interval Reachability Analysis: 
% Bounding Trajectories of Uncertain Systems with Boxes for Control and 
% Verification". Springer, 2021. DOI: 10.1007/978-3-030-65110-7

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

addpath('Utilities');

%% Logger configuration
% The logger is used to output information at various criticality levels.

% Configure logging level.
% Setting to a particular logging level will allow messages
% for equal or higher levels to be printed to the console.
% For example, setting parameters.log_level = logger.logger.WARNING
% allows warning and error messages to be printed.
%   ALL     = 0
%   RUNTIME = 1
%   INFO    = 2
%   WARNING = 3
%   ERROR   = 4
%   OFF     = 5
parameters.log_level = logger.logger.ALL;

% Log to command window or to file.
%   NaN: to command window
%   file name: will append to that file.
%       e.g. parameters.log_location = "tira.log"
parameters.log_location = NaN;

% *The 3 commands below do not need to be modified by the user*
[log, delete_logger] = logger.get_logger();
log.set_log_level(parameters.log_level);
log.set_log_location(parameters.log_location);

%% ODE configuration
% Allows the user to select (amongst built-in solvers) or define their own
% ODE solver to be used in the methods that need one.

% Select or define ODE solver.
% Must have the same function signature as built-in ODE solvers invoked as:
%   parameters.ode_solver(odefun, tspan, y0, parameters.ode_options)
% This format excludes the fully implicit solvers ode15i and decic.
% Inputs:
%   odefun: function which takes time and state and returns derivatives.
%       dydt = odefun(t, y)
%   tspan: time interval of the integration. E.g. [t_init, t_final].
%   y0: initial conditions.
%   options: option structure.
% Outputs:
%   t: evaluation points.
%   y: solutions at corresponding time points in t.
parameters.ode_solver = @ode45;

% Options structure that will be passed to calls of ode_solver.
% If using a Matlab ODE solver, see 
% https://www.mathworks.com/help/matlab/math/summary-of-ode-options.html
% for available options and how to construct the options variable.
parameters.ode_options = {};

%% Choice of the over-approximation method

% 1 - Continuous-time monotonicity
%     This choice will check if the system is monotone, then apply an
%     over-approximation method specific to monotone systems.
%     Pick option 3 for a slightly slower version but which does not
%     require checking the monotonicity satisfaction, yet will still
%     provide the same results if the system is indeed monotone.
% *Requirements:
% - Continuous-time system
% - One of the following 2 items:
%     Signs of the Jacobian matrices in ./Input_files/UP_Jacobian_Signs.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m

% 2 - Continuous-time contraction and growth bound
%     This method relies on using or defining a function that bounds the
%     growth or contraction on each dimension of the system.
% *Requirements:
% - Continuous-time system
% - One of the following 3 items:
%     Growth bound function in ./Input_files/UP_Growth_Bound_Function.m
%     Contraction matrix or scalar in ./Input_files/UP_Contraction_Matrix.m
%     Bounds of the state Jacobian in ./Input_files/UP_Jacobian_Bounds.m
%     *The last 2 cases also require dynamics with additive input:
%     System_description(t,x,p)=System_description(t,x,0)+p
%     which is automatically checked in ./Utilities/Growth_bound_choice.m

% 3 - Continuous-time mixed-monotonicity
%     This method is an extension of method 1 to any system with bounded
%     Jacobian matrices.
% *Requirements:
% - Continuous-time system
% - One of the following 2 items:
%     Signs of the Jacobian matrices in ./Input_files/UP_Jacobian_Signs.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m

% 4 - Continuous-time sampled-data mixed-monotonicity
%     This method applies the discrete-time mixed-monotonicity method (5)
%     to the sampled-data version of a continuous-time system.
%     The requirements on the Jacobian of the discrete-time system thus
%     translate into requirements on the sensitivity matrices of the 
%     continuous-time system.
% *Requirements:
% - Continuous-time system
% - Input p constant over the considered time range [t_init,t_final].
% - Optional user-provided inputs amongst the following items:
%     Signs of the sensitivity matrices in ./Input_files/UP_Sensitivity_Signs.m
%     Bounds of the sensitivity matrices in ./Input_files/UP_Sensitivity_Bounds.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m
%     Bounds of the 2nd order Jacobian matrices in ./Input_files/UP_Jacobian_2ndOrder_Bounds.m
%     Jacobian functions in ./Input_files/UP_Jacobian_Function.m
%     *If none of these optional input is provided, the method can still
%     be called using default submethod 4-0 "Sampling and falsification".

% 5 - Discrete-time mixed-monotonicity
%     Similar to method 3 but for discrete-time systems.
% *Requirements:
% - Discrete-time system
% - One of the following 2 items:
%     Signs of the Jacobian matrices in ./Input_files/UP_Jacobian_Signs.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m

% 6 - Quasi-Monte Carlo
%     This method simulates the system over a grid of points and combines
%     it with an error-bound function to over-approximate the reachable set.
% *Requirements for discrete-time systems:
% - Error-bound function in ./Input_files/UP_Error_Bound_Function.m
% *Requirements for continuous-time systems:
% - One of the following 4 items:
%     Error-bound function in ./Input_files/UP_Error_Bound_Function.m
%     Lipschitz constant in ./Input_files/UP_Lipschitz_Constant.m
%     Contraction matrix or scalar in ./Input_files/UP_Contraction_Matrix.m
%     Bounds of the Jacobian matrices in ./Input_files/UP_Jacobian_Bounds.m
%     *The last 2 cases also require dynamics with additive input:
%     System_description(t,x,p)=System_description(t,x,0)+p
%     which is automatically checked in ./Utilities/Error_bound_choice.m

% 7 - Monte Carlo
%     This method simulates the system over uniformly sampled points to 
%     generate an approximation interval with a desired probabilistic 
%     accuracy and confidence.
%     *The resulting interval approximation is not guaranteed to be 
%     an OVER-approximation of the reachable set.
% *This method has no requirement and is applicable to both
% discrete-time and continuous-time systems.

% NaN - (no specific method requested)
%     This lets the function TIRA.m automatically select the best available 
%     method based on system properties and provided additional system 
%     description in the files of the './Input_files' folder

parameters.OA_method = NaN;

    %% (2): Contraction/Growth bound submethod choice

    % How is the growth bound function obtained?
    %   2-1: Read it from user-provided file 
    %       Requires Input_files/UP_Growth_Bound_Function.m
    %   2-2: Create it from user-provided contraction matrix or scalar
    %       Requires Input_files/UP_Contraction_Matrix.m
    %   2-3: Create it from user-provided Jacobian bounds 
    %       Requires Input_files/UP_Jacobian_Bounds.m)
    %   NaN: let the function Utilities/Growth_bound_choice.m select the 
    %       best available method (in the above order) based on provided 
    %       files in './Input_files' folder
    parameters.growth_bound_method = NaN;

    %% (4): Sampled-data mixed-monotonicity: submethod choice

    % How are bounds on the sensitivity matrices obtained?
    %	4-1: Read sensitivity bounds from user-provided files 
    %       Requires either Input_files/UP_Sensitivity_Bounds.m 
    %       or Input_files/UP_Sensitivity_Signs.m
    %	4-2: Interval Arithmetics 
    %       Requires Input_files/UP_Jacobian_Function.m
    %   4-3: 2nd order Sensitivity 
    %       Requires Input_files/UP_Jacobian_Bounds.m 
    %       and Input_files/UP_Jacobian_2ndOrder_Bounds.m
    %	4-0: Sampling and falsification 
    %       No additional requirement.
    %   NaN: let the function Utilities/Sensitivity_Bounds.m select the 
    %       best available method based on provided files in 
    %       './Input_files' folder
    %       (Default if no file is provided: 4-0) 
    parameters.sensitivity_bounds_method = NaN;

        %% (4-2) and (4-3): Parameters for Taylor series of interval matrix exponential
        
        % Initial order for the truncation of the Taylor series 
        %   (may be interatively increased if results are not satisfactory)
        parameters.taylor_order = 100;
        
        % Stopping condition: 
        %   max relative error allowed between two iterations 
        %   (to do a single iteration, set: taylor_tolerance = NaN)
        parameters.taylor_tolerance = 1e-3;
        
        % Taylor order increase at each iteration
        parameters.taylor_increment = 100;
        
        % Scaling factor used in Scaling and Squaring method 
        %   (set to NaN to let the code pick it automatically) 
        parameters.scaling_factor = NaN;
               
        %% (4-3): Sampling parameters for the 2nd Order Sensitivity method
        
        % Number of sample points for each state dimension
        %   (needs to be a vector of same length as the state variable)
        parameters.sensitivity_samples_per_state_dimension = 2 * ones(n_x, 1);
        
        % Number of sample points for each input dimension
        %   (needs to be a vector of same length as the input variable)
        parameters.sensitivity_samples_per_input_dimension = 2 * ones(n_p, 1);
        
        %% (4-3) and (4-0): Solvers of the sensitivity equations
        %   0: Euler method on the sensitivity definition
        %   1: Euler method on Jacobian + ODE solver on sensitivity ODE
        %   2: ODE solver on sensitivity ODE using explicit Jacobian function
        %       Requires Input_files/UP_Jacobian_Function.m
        %   NaN: let the function Utilities/Sensitivity_Bounds.m select the
        %       best available solver based on provided files in the 
        %       './Input_files' folder
        %       (Default if no file is provided: 0) 
        parameters.sensitivity_equation_solver = NaN;
        
        % Parameters for Euler approximations in Sensitivity_solver.m
        %   when sensitivity_equation_solver is chosen as 0 or 1
        parameters.euler_epsilon = 1e-4;         % Initial step size
        parameters.euler_tolerance = 1e-10;      % Error tolerance between iterations (stopping condition of the loop) 
        parameters.euler_max_iter = 4;           % Max number of iterations (each iteration divides the step by 2)

        %% (4-0): Sampling and falsification parameters
        
        % Sample pairs in [x_low,x_up]*[p_low,p_up] for which sensitivities
        %   are evaluated in function Utilities/Sensitivity_sampling.m
        %   NaN: let the sampling function pick a default value
        %       (sensitivity_sampling_per_dimension^(n_x+n_p))
        %   1: Skip and only do falsification. 
        parameters.sensitivity_sampling_pairs = NaN;
        
        % Number of sample pairs per state and input dimension to be used
        %   in Utilities/Sensitivity_sampling.m 
        %   when 'sensitivity_sampling_pairs = NaN'
        parameters.sensitivity_sampling_per_dimension = 2;

        % Max number of falsification iterations in Utilities/Sensitivity_falsification.m
        %   Inf: infinite
        %   0: Skip falsification
        parameters.falsification_max_iter = Inf;     

    %% (6): Quasi-Monte Carlo, sampling parameters and submethod choice
    
    % Number of sample points for each state dimension
    %   (needs to be a vector of same length as the state variable)
    parameters.QMC_samples_per_state_dimension = 2 * ones(n_x, 1);

    % Number of sample points for each input dimension
    %   (needs to be a vector of same length as the input variable)
    parameters.QMC_samples_per_input_dimension = 2 * ones(n_p, 1);
    
    % How is the error bound function obtained?
    %   6-1: Read it from user-provided file 
    %       Requires Input_files/UP_Error_Bound_Function.m
    %   6-2: Create it from user-provided contraction matrix or scalar
    %       Requires Input_files/UP_Contraction_Matrix.m
    %   6-3: Create it from user-provided Jacobian bounds 
    %       Requires Input_files/UP_Jacobian_Bounds.m)
    %   6-4: Create it from user-provided Lipschitz constant
    %       Requires Input_files/UP_Lipschitz_Constant.m
    %   NaN: let the function Utilities/Error_bound_choice.m select the 
    %       best available method (in the above order) based on provided 
    %       files in './Input_files' folder
    parameters.error_bound_method = NaN;
    
    %% (7): Monte Carlo
    
    % Accuracy of the resulting interval approximation
    %   (define epsilon in (0, 1) to obtain a (1-epsilon)-accuracy)
    parameters.epsilon = 0.05;
    
    % Confidence in obtaining a (1-epsilon)-accuracy
    %   (define delta in (0, 1) to obtain a (1-delta)-confidence)
    parameters.delta = 1e-9;  

%% Default values
% parameters.log_level = logger.logger.ALL;
% parameters.log_location = NaN;
% parameters.ode_solver = @ode45;
% parameters.ode_options = {};

% parameters.OA_method = NaN;
% parameters.growth_bound_method = NaN;
% parameters.sensitivity_bounds_method = NaN;
% parameters.taylor_order = 100;      
% parameters.taylor_tolerance = 1e-3; 
% parameters.taylor_increment = 100;
% parameters.scaling_factor = NaN;
% parameters.sensitivity_samples_per_state_dimension = 2 * ones(n_x, 1);
% parameters.sensitivity_samples_per_input_dimension = 2 * ones(n_p, 1);
% parameters.sensitivity_equation_solver = NaN;
% parameters.euler_epsilon = 1e-4;         
% parameters.euler_tolerance = 1e-10;      
% parameters.euler_max_iter = 4;
% parameters.sensitivity_sampling_pairs = NaN;            
% parameters.sensitivity_sampling_per_dimension = 2;      
% parameters.falsification_max_iter = Inf;   
% parameters.QMC_samples_per_state_dimension = 2 * ones(n_x, 1);
% parameters.QMC_samples_per_input_dimension = 2 * ones(n_p, 1);
% parameters.error_bound_method = NaN;
% parameters.epsilon = 0.05; % 95% accuracy
% parameters.delta = 1e-9; % (1-10^(-9)) confidence in accuracy (1-epsilon)
