%% User-provided growth bound function handle for a continuous-time system
% Used in:
% - Utilities/Growth_bound_choice.m
%       for the over-approximation method 2 in
%       OA_methods/OA_2_CT_Contraction_growth_Bound.m

% The main requirements on its definition for the over-approximation to be 
% applicable are provided at the end of this file, or in more details in 
% the paper below.

% Source paper:
% G. Reissig, A. Weber and M. Rungger, "Feedback refinement relations for 
% the synthesis of symbolic controllers". IEEE Transactions on Automatic 
% Control v. 62(4), pp. 1781-1796, 2017. DOI: 10.1109/TAC.2016.2593947

% List of inputs
%   t_init: initial time
%   t_final: time at which the reachable set is approximated
%   [x_low,x_up]: interval of initial states (at time t_init)
%   [p_low,p_up]: interval of allowed input values

% List of outputs
%   GB_handle: function handle of the growth bound
%       inputs: positive time, positive state vector, positive input vector
%       output: positive state vector (growth or contraction of the system)

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function GB_handle = UP_Growth_Bound_Function(t_init,t_final,x_low,x_up,p_low,p_up)
n_x = length(x_low);
n_p = length(p_low);

%% Default values as NaN (not a number)
GB_handle = @(t,x,p) NaN(n_x,1);

%% User-provided growth bound handle

% Requirements on the definition of GB_handle
% (using componentwise inequalities)
%   1) GB_handle is a function handle from R+*R+^n_x*R+^n_p to R+^n_x
% 
%   2) x>=x', p>=p' => GB_handle(t,x,p)>=GB_handle(t,x',p')
% 
%   3) Let x(t_final;t_init,x0,p) be the solution of System_description(t,x,p)
%   at time t_final, starting from x0 at t_init and with constant input p.
%   Let x0_c and p_c be the centers of [x_low,x_up] and [p_low,p_up].
%   Then we need, for all x0 in [x_low,x_up] and p in [p_low,p_up]:
%   abs(x(t_final;t_init,x0,p)-x(t_final;t_init,x0_c,p_c))
%       <= GB_handle(t_final-t_init,abs(x0-x0_c),abs(p-p_c))


% GB_handle = @(t,x,p) ...
global system_choice

switch system_choice
    case 99
        %% Unicycle with additive disturbance on (x,y)  -- Eq. (11)
        % State x = [x; y; theta; v]
        % Input p = [u1; u2; d1; d2], ZOH over [t_init, t_final]
        %
        % We mirror the structure of case 11:
        %   - build a contraction matrix C from Jacobian bounds wrt x
        %   - build an input-bound vector p_tilde that upper-bounds input effects
        %
        % IMPORTANT: Our input intervals are symmetric:
        %   u1 in [-1,1], u2 in [-2,2], d1,d2 in [-0.3,0.3]
        %   => input_center = 0, so we can use the radii directly.

        % center of the input box (should be zero with symmetric bounds)
        input_center = (p_low + p_up)/2;

        % Jacobian bounds wrt x over [x_low,x_up] at input_center
        [J_x_low, J_x_up, ~, ~] = UP_Jacobian_Bounds( ...
            t_init, t_final, x_low, x_up, input_center, input_center);

        % Contraction matrix C (same pattern as case 11)
        C = max(abs(J_x_low), abs(J_x_up));   % off-diagonal: abs upper bound
        C(1:n_x+1:end) = diag(J_x_up);        % diagonal: take upper bound directly

        % ---- Input-effect bound vector p_tilde (maps p to state rates) ----
        % Dynamics: 
        %   xdot = v*cos(theta) + d1
        %   ydot = v*sin(theta) + d2
        %   thetadot = u1
        %   vdot = u2
        %
        % With symmetric input box, the radius of each input is:
        p_rad = abs(p_up - input_center);   % = [1; 2; 0.3; 0.3] for your case
        % Map (u1,u2,d1,d2) -> (xdot,ydot,thetadot,vdot):
        %   d1 affects xdot, d2 affects ydot, u1 affects thetadot, u2 affects vdot.
        p_tilde = [ p_rad(3);     % |d1| bound -> xdot
                    p_rad(4);     % |d2| bound -> ydot
                    p_rad(1);     % |u1| bound -> thetadot
                    p_rad(2) ];   % |u2| bound -> vdot

        % Growth bound handle (identical structure to case 11)
        GB_handle = @(t, x, p) expm(C*t)*x + ...
            integral(@(s) expm(C*s)*p_tilde, 0, t, 'ArrayValued', true);
end


