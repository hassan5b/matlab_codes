%% User-provided Lipschitz constant for a continuous-time system
% Used in:
% - Utilities/Error_bound_choice.m
%       for the definition of an error bound function to be used in
%       over-approximation method 6 in
%       OA_methods/OA_6_DT_CT_Quasi_Monte_Carlo.m

% List of inputs:
%   t_init: initial time
%   t_final: time at which the reachable set is approximated
%   [x_low, x_up]: interval of initial states (at time t_init).
%   [p_low, p_up]: interval of allowed input values.

% List of outputs:
%   L: Lipschitz constant for the system

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function L = UP_Lipschitz_Constant(t_init, t_final, x_low, x_up, p_low, p_up)
n_x = length(x_low);
n_p = length(p_low);

%% Default values as NaN (not a number)
L = NaN;

%% User-provided Lipschitz constant L

% Requirements on the definition of L
%   For the continuous-time system whose vector field is defined as
%   System_description(t,x,p), we need for all states x,y and inputs p,q
%       norm(System_description(t,x,p)-System_description(t,y,q),Inf) 
%           <= L*norm([x-y;p-q],Inf)

global system_choice;

switch system_choice
    case 6
        L = 3;
end
