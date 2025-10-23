%% Main test file calling most over-approximation methods and submethods
% to ensure that future changes or additions to the toolbox
% do not cause unexpected damage. The tests currently implemented are:
% - Calling each over-approximation method 
% - Checking that the correct method is selected when none is specified
% - Checking that the correct submethod is selected when none is specified

% This file writes down the success or failure of each test in the file
% tira_test_output.log

% *Note: the failure or success of each test is determined by scanning for
% the expected logger outputs (log.info('...')) of each method.
% Modifications in the content of these logging strings in TIRA will need 
% to be reflected here as well.

% Authors:  
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

%% Initialization
close all;
clear;
assert(isfile('Run_tests.m'),'Run_tests.m needs to be run from the Test_suite folder')
addpath('../OA_methods');
addpath('../Utilities');

% Create temporary file for TIRA logging.
log_name = 'tira_test.log';
fopen(log_name, 'w'); % Clear log. Must be cleared.

global system_choice;
global u;

global g_sensitivity_bounds_method;
global g_error_bound_method;
g_sensitivity_bounds_method = NaN;
g_error_bound_method = NaN;

% fopen("tira_test_output.log", 'w'); % Clear output. Do not clear if old test results need to be kept.
test_output_log = logger.logger();
test_output_log.set_log_level(logger.logger.ALL);
test_output_log.set_log_location("tira_test_output.log");
test_output_log.info('Beginning new log:\n\n');

%% Verify directly calling an over-approximation method calls the desired method for valid input

test_output_log.info('Test directly calling over-approximation methods:\n');

% Over-approximation method 1: Continuous-time monotonicity
% Protein interaction (continuous-time)
system_choice = 8;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 1);
check_log_contains(test_output_log, log_name, 'continuous-time monotonicity');
        
% Over-approximation method 2: Continuous-time componentwise contraction and growth bound
% Unicycle (continuous-time)
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final, u] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 2);
check_log_contains(test_output_log, log_name, 'contraction/growth bound');

% Over-approximation method 3: Continuous-time mixed-monotonicity
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 3);
check_log_contains(test_output_log, log_name, 'continuous-time mixed-monotonicity');

% Over-approximation method 4: Continuous-time sampled-data mixed-monotonicity
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 4);
check_log_contains(test_output_log, log_name, 'sampled-data mixed-monotonicity');

% Over-approximation method 5: Discrete-time mixed-monotonicity
% Traffic diverge 3 link (discrete-time)
system_choice = 4;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA(t_init, x_low, x_up, p_low, p_up, 5);
check_log_contains(test_output_log, log_name, 'discrete-time mixed-monotonicity');

% Over-approximation method 6: Discrete or Continuous-time quasi-Monte Carlo
% Unicycle (continuous-time)
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 6);
check_log_contains(test_output_log, log_name, 'quasi-Monte Carlo');

% Simple unstable model (discrete-time)
system_choice = 14;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA(t_init, x_low, x_up, p_low, p_up, 6);
check_log_contains(test_output_log, log_name, 'quasi-Monte Carlo');

% Over-approximation method 7: Discrete or Continuous-time Monte Carlo
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 7);
check_log_contains(test_output_log, log_name, 'Monte Carlo');

% Traffic diverge 3 link (discrete-time)
system_choice = 4;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA(t_init, x_low, x_up, p_low, p_up, 7);
check_log_contains(test_output_log, log_name, 'Monte Carlo');

%% Check that the correct over-approximation method is selected when none is specified
test_output_log.info('Check that the correct over-approximation method is selected when none is specified:\n');

% Unicycle (continuous-time)
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up);
check_log_contains(test_output_log, log_name, 'contraction/growth bound');

% Traffic diverge 3 link (discrete-time)
system_choice = 4;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA(t_init, x_low, x_up, p_low, p_up);
check_log_contains(test_output_log, log_name, 'discrete-time mixed-monotonicity');

% Protein interaction (continuous-time)
system_choice = 8;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up);
check_log_contains(test_output_log, log_name, 'continuous-time mixed-monotonicity');

%% Check that the correct submethod is selected when none is specified
test_output_log.info('Check that the correct submethod is selected when none is specified\n');

% Continuous-time componentwise contraction and growth bound
% Defined growth bound function
system_choice = 11;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 2);
check_log_contains(test_output_log, log_name,...
    'contraction/growth bound',...
    'Extract growth bound function handle from Input_files/UP_Growth_Bound_Function.m');

% Contraction matrix
system_choice = 6;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 2);
check_log_contains(test_output_log, log_name,...
    'contraction/growth bound',...
    'Extract growth bound function handle from contraction matrix in Input_files/UP_Contraction_Matrix.m');

% Jacobian bounds
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 2);
check_log_contains(test_output_log, log_name,...
    'contraction/growth bound',...
    'Extract growth bound function handle from state Jacobian bounds in Input_files/UP_Jacobian_Bounds.m');

% Sampled-data mixed-monotonicity
% Interval arithmetics
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 4);
check_log_contains(test_output_log, log_name,...
    'sampled-data mixed-monotonicity',...
    'Extract sensitivity bounds with interval arithmetics');

% Second-order sensitivities
system_choice = 12;
g_sensitivity_bounds_method = 3;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 4);
check_log_contains(test_output_log, log_name,...
    'sampled-data mixed-monotonicity',...
    'Extract sensitivity bounds with second-order sensitivity method.');

% Sampling and falsification
system_choice = 6;
g_sensitivity_bounds_method = NaN;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 4);
check_log_contains(test_output_log, log_name,...
    'sampled-data mixed-monotonicity',...
    'Extract sensitivity bounds with the sampling and falsification method.');

% Quasi-Monte Carlo
% Error-bound function
system_choice = 6;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 6);
check_log_contains(test_output_log, log_name,...
    'quasi-Monte Carlo',...
    'Extract error bound function handle from Input_files/UP_Error_Bound_Function.m');

% Contraction matrix
system_choice = 7;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 6);
check_log_contains(test_output_log, log_name,...
    'quasi-Monte Carlo',...
    'Extract error bound function from contraction matrix');

% Jacobian bounds
system_choice = 1;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 6);
check_log_contains(test_output_log, log_name,...
    'quasi-Monte Carlo',...
    'Extract error bound function from state Jacobian bounds in Input_files/UP_Jacobian_Bounds.m');

% Lipschitz constant
system_choice = 6;
g_error_bound_method = 4;
[x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final] = setup_system(system_choice);
TIRA([t_init, t_final], x_low, x_up, p_low, p_up, 6);
check_log_contains(test_output_log, log_name,...
    'quasi-Monte Carlo',...
    'Extract error bound function from Lipschitz constant');

%% Cleanup
% Delete temporary file.
delete(log_name);

%% Scan the logger output of each method called and return "success" or "failure" log line
% This is highly dependent on the format of the log strings currently used
% in TIRA.m and the Utilities folder. Changes in the strings of these
% log.info command need to be reflected here as well.
function check_log_contains(logger, log_name, string, string2)
    % Read the logger output for the method that was called
    text = fileread(log_name);
    
    % Look in the logger output for a specific string "Call of the <method-name>"
    %   that is only displayed when the method could be successfully called
    % Write in the test-suite logger success or failure depending on
    %   whether this string could be found
    if ~contains(text, ['Call of the ', string])
        % Continue to next test if errors.
        try
            logger.error(['Failed to call ', string, '\n']);
        catch
        end
    else
        logger.info(['    Succeeded in calling ', string, '\n']);
    end
    
    % When a second strong is provided (for a submethod name), do the same
    % check to see if this submethod was sucessfully called
    if nargin == 4
        if ~contains(text, string2)
            % Continue to next test if errors.
            try
                logger.error(['Failed to ', string2, '\n']);
            catch
            end
        else
            logger.info(['    Succeeded in ', string2, '\n']);
        end
    end
    
    % Clear the contents of the log file for this method
    %   (this doesn't affect the main logger of the test suite)
    fopen(log_name, 'w');
end

%% Read the pre-defined reachability problem definition for the requested system index
function [x_low, x_up, n_x, p_low, p_up, n_p, t_init, t_final, u] = setup_system(system_id)
    x_low = NaN;
    x_up = NaN;
    n_x = NaN;
    p_low = NaN;
    p_up = NaN;
    n_p = NaN;
    t_init = NaN;
    t_final = NaN;
    u = NaN;
    switch system_id
        case 1
            x_up = [1;1;pi/8];
            x_low = [0;0;-pi/8];
            n_x = length(x_up);
            u = [0.25;0.15];
            p_up = [0.05;0.05;0.03];
            p_low = -p_up;
            n_p = length(p_up);
            t_init = 0;
            t_final = 4;
        case 4
            x_low = [150;180;100];
            x_up = [200;300;220];
            n_x = length(x_low);
            p_low = 40;
            p_up = 60;
            n_p  = length(p_low);
            t_init = 0;
        case 6
            %% Time-varying system (continuous-time)
            % Interval of initial states (defined by two column vectors)
            x_low = [1;2;3];
            x_up = [2;4;4];
            n_x = length(x_low);

            % Interval of allowed input values (defined by two column vectors)
            p_up = 4*ones(n_x,1);
            p_low = -p_up;
            n_p  = length(p_low);

            % Time interval
            t_init = 1;
            t_final = 2;
        case 7
            %% Bicycle (continuous-time)
            % Interval of initial states (defined by two column vectors)
            x_up = [1;1;pi/8];
            x_low = [0;0;-pi/8];
            n_x = length(x_up);

            % Control
            u = [0.25;0.15];

            % Interval of allowed input values (disturbances)
            p_up = 0.3*ones(n_x,1);
            p_low = -p_up;
            n_p  = length(p_low);

            % Time interval
            t_init = 0;
            t_final = 1;
        case 8
            x_low = [0;0];
            x_up = [1;1];
            n_x = length(x_low);
            p_up = 0;
            p_low = 0;
            n_p  = length(p_low);
            t_init = 0;
            t_final = 1;
        case 11
            %% Pursuer-evader game with 2 Dubin's vehicles (continuous-time)
            % Interval of initial states (defined by two column vectors)
            x_up = [1;1;pi/8;2.5;2.5;3*pi/8];
            x_low = [0;0;-pi/8;2;2;pi/4];
            n_x = length(x_up);

            % Control
            u = [0.25;0.15];

            % Interval of allowed input values (disturbances)
            p_up = [0.5;0.3];
            p_low = -p_up;
            n_p = length(p_up);

            % Time interval
            t_init = 0;         % Initial time
            t_final = 1;        % Final time
        case 12
            %% Simple oscillator model (continuous-time)
            % Interval of initial states (defined by two column vectors)
            x_low = [0;0];
            x_up = [1;1];
            n_x = length(x_low);

            % Interval of allowed input values
            % (system has no input: set both p_up and p_low to 0)
            p_up = [0;0];
            p_low = [0;0];
            n_p  = length(p_low);

            % Time interval
            t_init = 0;
            t_final = 1;
        case 14
            %% Simple unstable model (discrete-time)
            % Interval of initial states (defined by two column vectors)
            x_low = [-1; -1];
            x_up = [1; 1];
            n_x = length(x_low);
            
            % Interval of allowed input values
            p_up = [1; 1];
            p_low = [-1; -1];
            n_p = length(p_low);
            
            % Time interval
            t_init = 0;            
            
    end
end