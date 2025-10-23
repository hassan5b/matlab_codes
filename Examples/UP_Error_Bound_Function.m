%% User-provided error bound function handle for continuous-time or discrete-time systems
% Used in:
% - Utilities/Error_bound_choice.m
%       for the over-approximation method 6 in
%       OA_methods/OA_6_DT_CT_Quasi_Monte_Carlo.m

% Source paper:
% P.-J. Meyer, A. Devonport and M. Arcak, "Sampling-Based Methods", 
% In: "Interval Reachability Analysis: Bounding Trajectories of Uncertain 
% Systems with Boxes for Control and Verification". Springer, 2021
% DOI: 10.1007/978-3-030-65110-7_7

% List of inputs:
%   time_vector: time information for the reachability analysis
%       1D vector for discrete-time system: time_vector = t_init
%       2D vector for continuous-time system: time_vector=[t_init, t_final]
%       t_init: initial time
%       t_final: if present, time at which  reachable set is approximated.
%           For discrete-time systems, t_final will be one step from t_init
%   [x_low, x_up]: interval of initial states (at time t_init).
%   [p_low, p_up]: interval of allowed input values.

% List of outputs:
%   error_bound_handle: function handle of error bound.
%       inputs: time, infinity norm of error between state, and infinity
%           norm of error between input.
%       output: Bound on infinity norm error between evolved states
%           given input errors.

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function error_bound_handle = UP_Error_Bound_Function(time_vector,x_low,x_up,p_low,p_up)
n_x = length(x_low);
n_p = length(p_low);

%% Default values as NaN (not a number)
error_bound_handle = @(time_vector, state_error, input_error) NaN;

%% User-provided error bound handle

% Requirements on the definition of error_bound_handle
% (using componentwise inequalities)
%   1) error_bound_handle is a function handle 
%       from R^2*R+*R+ to R+ for continuous-time systems
%       from R*R+*R+ to R+ for discrete-time systems
% 
%   2) x>=x', p>=p' => error_bound_handle(t,x,p)>=error_bound_handle(t,x',p')
% 
%   3) For a continuous-time system, let x(t_final;t_init,x0,p) be the 
%       solution of System_description(t,x,p) at time t_final, starting 
%       from x0 at t_init and with constant input p.
%       Then we need, for all x,y in [x_low,x_up] and p,q in [p_low,p_up]:
%       norm(x(t_final;t_init,x,p)-x(t_final;t_init,y,q),Inf)
%           <= error_bound_handle([t_init,t_final],norm(x-y,Inf),norm(p-q,Inf))
%
%   4) For a discrete-time system, let succ(t_init,x0,p) be the one-step
%       successor of system System_description(t,x,p) with initial state x0
%       and input p.
%       Then we need, for all x,y in [x_low,x_up] and p,q in [p_low,p_up]:
%       norm(succ(t_init,x,p)-succ(t_init,y,q),Inf)
%           <= error_bound_handle(t_init,norm(x-y,Inf),norm(p-q,Inf))

global system_choice

switch system_choice
    case 6
        %% Time-varying system (continuous-time)
        % This is a simple linear periodic time-varying system, 
        % without a particular physical interpretation attached to it.
        % dx = A(t)*x+p
        error_bound_handle = @(time_vector, state_error, input_error)...
                exp(3 * (time_vector(2) - time_vector(1))) * (state_error + input_error);
end

