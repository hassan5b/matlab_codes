%% User-provided bounds on the Jacobian matrices for continuous-time or discrete-time systems
% Used in:
% - OA_methods/TIRA.m
%       to check the requirements for over-approximation methods 1,3,5 in
%       OA_methods/OA_1_CT_Monotonicity.m
%       OA_methods/OA_3_CT_Mixed_Monotonicity.m
%       OA_methods/OA_5_DT_Mixed_Monotonicity.m
% - Utilities/Growth_bound_choice.m
%       for the definition of a growth bound function to be used in
%       over-approximation method 2 in
%       OA_methods/OA_2_CT_Contraction_growth_Bound.m
% - Utilities/Sensitivity_bounds_choice.m
%       to compute bounds on the sensitivity matrices of a continuous-time
%       system, within submethods (4-2) and (4-3), to be used in 
%       over-approximation method 4 in
%       OA_methods/OA_4_CT_Sampled_data_MM.m
% - Utilities/Error_bound_choice.m
%       for the definition of an error bound function to be used in
%       over-approximation method 6 in
%       OA_methods/OA_6_DT_CT_Quasi_Monte_Carlo.m

% The user can either provider global bounds, or a function of the input
% arguments (t_init,t_final,x_low,x_up,p_low,p_up) returning local bounds

% Jacobian definitions:
%   to states:  J_x(t) = d(System_description(t,x,p))/dx
%   to inputs:  J_p(t) = d(System_description(t,x,p))/dp
% Note that Jacobians may contain elements inf or -inf, but not all methods
% may support this.

% List of inputs
%   t_init: initial time
%   t_final: time at which the reachable set is approximated (for continuous-time system only)
%       for discrete-time system, a dummy value can be provided
%   [x_low,x_up]: interval of initial states (at time t_init)
%   [p_low,p_up]: interval of allowed input values

% List of outputs
%   [J_x_low,J_x_up]: bounds of the Jacobian with respect to the state
%   [J_p_low,J_p_up]: bounds of the Jacobian with respect to the input

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function [J_x_low,J_x_up,J_p_low,J_p_up] = UP_Jacobian_Bounds(t_init,t_final,x_low,x_up,p_low,p_up)
n_x = length(x_low);
n_p  = length(p_low);

%% Default values as NaN (not a number)
J_x_low = NaN(n_x);
J_x_up = NaN(n_x);
J_p_low = NaN(n_x,n_p);
J_p_up = NaN(n_x,n_p);

%% User-provided Jacobian bounds
% Can be either global bounds
% or local bounds depending on inputs: t_init,t_final,x_low,x_up,p_low,p_up

% If System_description.m has no input variable 'p', uncomment below:
% J_p_low = zeros(n_x,n_p);
% J_p_up = zeros(n_x,n_p);
disp(['DBG UB file = ', mfilename('fullpath')])
global system_choice
disp(['DBG system_choice in UB = ', num2str(system_choice)])


switch system_choice
   case 99
    % Unicycle + disturbance on (x,y): x=[x;y;theta;v], p=[u1;u2;d1;d2]
    % Continuous-time, ZOH over [t_init,t_final]
    disp('DBG: ENTERED UP_Jacobian_Bounds case 99')

    T = t_final - t_init;
    n_x = length(x_low); n_p = length(p_low);

    % --- init outputs with FINITE numbers (no NaNs) ---
    J_x_low = zeros(n_x);  J_x_up  = zeros(n_x);
    J_p_low = zeros(n_x,n_p); J_p_up = zeros(n_x,n_p);

    % --- bound theta, v over the step (ZOH on u1,u2) ---
    thL = x_low(3) + min(0, T * p_low(1));   % theta^- = theta0 + T * u1^-
    thU = x_up(3)  + max(0, T * p_up(1));    % theta^+ = theta0 + T * u1^+
    vL  = x_low(4) + min(0, T * p_low(2));   % v^-     = v0 + T * u2^-
    vU  = x_up(4)  + max(0, T * p_up(2));    % v^+     = v0 + T * u2^+
    vabs = max(abs([vL, vU]));

    % Extremely safe trig bounds (always finite)
    [cmn,cmx] = cos_bounds(thL, thU);
    [smn,smx] = sin_bounds(thL, thU);
    
    J_x_low(1,3) = -vabs * smx;  J_x_up(1,3) = -vabs * smn;
    J_x_low(2,3) =  vabs * cmn;  J_x_up(2,3) =  vabs * cmx;
    
    J_x_low(1,4) =  cmn;         J_x_up(1,4) =  cmx;
    J_x_low(2,4) =  smn;         J_x_up(2,4) =  smx;


    % --- Jacobian wrt inputs p = [u1;u2;d1;d2] (constant mapping) ---
    % d1->xdot, d2->ydot, u1->thetadot, u2->vdot
    J_p_low(1,3) = 1;  J_p_up(1,3) = 1;    % ∂xdot/∂d1
    J_p_low(2,4) = 1;  J_p_up(2,4) = 1;    % ∂ydot/∂d2
    J_p_low(3,1) = 1;  J_p_up(3,1) = 1;    % ∂thetadot/∂u1
    J_p_low(4,2) = 1;  J_p_up(4,2) = 1;    % ∂vdot/∂u2

    % (all other entries remain zero and finite)

end

% tight cos/sin bounds on [a,b]
function [cmin,cmax] = cos_bounds(a,b)
    if b < a, [a,b] = deal(b,a); end
    if b-a >= 2*pi, cmin=-1; cmax=1; return; end
    cmin = min(cos(a),cos(b)); cmax = max(cos(a),cos(b));
    % maxima at 2kπ, minima at π+2kπ
    if floor(b/(2*pi)) - ceil(a/(2*pi)) >= 0, cmax = 1; end
    if floor((b-pi)/(2*pi)) - ceil((a-pi)/(2*pi)) >= 0, cmin = -1; end
end
function [smin,smax] = sin_bounds(a,b)
    if b < a, [a,b] = deal(b,a); end
    if b-a >= 2*pi, smin=-1; smax=1; return; end
    smin = min(sin(a),sin(b)); smax = max(sin(a),sin(b));
    % maxima at π/2+2kπ, minima at -π/2+2kπ
    if floor((b - pi/2)/(2*pi)) - ceil((a - pi/2)/(2*pi)) >= 0, smax = 1; end
    if floor((b + pi/2)/(2*pi)) - ceil((a + pi/2)/(2*pi)) >= 0, smin = -1; end
end
