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

global system_choice
global u        % Some systems need an external control input

switch system_choice
    case 1
        %% Unicycle (continuous-time)
        time_step = t_final-t_init;
        n_x = length(x_low);
        J_x_low = zeros(n_x);
        J_x_up = zeros(n_x);
        J_p_low = eye(n_x);
        J_p_up = eye(n_x);

        % Possible orientation values over the next sampling period
        theta_low = x_low(3) + min(0,time_step*(u(2)+p_low(3)));
        theta_up = x_up(3) + max(0,time_step*(u(2)+p_up(3)));

        % Min/Max for cosine
        % (conditions only valid if x(3) initial set is in [-pi,pi])
        min_cos = -1;
        if ~((theta_low <= -pi && theta_up >= -pi) || (theta_low <= pi && theta_up >= pi))
            min_cos = min(cos(theta_low),cos(theta_up));
        end
        max_cos = 1;
        if ~((theta_low <= -2*pi && theta_up >= -2*pi) || (theta_low <= 0 && theta_up >= 0) || (theta_low <= 2*pi && theta_up >= 2*pi))
            max_cos = max(cos(theta_low),cos(theta_up));
        end
        
        % Min/Max for sine
        % (conditions only valid if x(3) initial set is in [-pi,pi])
        min_sin = -1;
        if ~((theta_low <= -5*pi/2 && theta_up >= -5*pi/2) || (theta_low <= -pi/2 && theta_up >= -pi/2) || (theta_low <= 3*pi/2 && theta_up >= 3*pi/2))
            min_sin = min(sin(theta_low),sin(theta_up));
        end
        max_sin = 1;
        if ~((theta_low <= -3*pi/2 && theta_up >= -3*pi/2) || (theta_low <= pi/2 && theta_up >= pi/2) || (theta_low <= 5*pi/2 && theta_up >= 5*pi/2))
            max_sin = max(sin(theta_low),sin(theta_up));
        end

        % Bounds of the partial derivatives
        if u(1) > 0
            J_x_low(1,3) = -u(1)*max_sin;
            J_x_up(1,3) = -u(1)*min_sin;
            J_x_low(2,3) = u(1)*min_cos;
            J_x_up(2,3) = u(1)*max_cos;
        else
            J_x_low(1,3) = -u(1)*min_sin;
            J_x_up(1,3) = -u(1)*max_sin;
            J_x_low(2,3) = u(1)*max_cos;
            J_x_up(2,3) = u(1)*min_cos; 
        end

    case 2
        %% Traffic diverge 3 links (continuous-time)
        % Defined as a system with additive input: dx=f(t,x)+p
        
        T = t_final-t_init;             % time period, in secondes
        v = 0.5;            % free-flow speed, in links/period
        w = 1/6;            % congestion-wave speed, in links/period

        J_x_low = 1/T*[-v 0 0;0 -(v+w) -w;0 -w -(v+w)];
        J_x_up = 1/T*[0 2*w 2*w;v/2 0 0;v/2 0 0];
        J_p_low = zeros(n_x);
        J_p_low(1,1) = 1;
        J_p_up = J_p_low;

    case 3
        %% Ship (continuous-time)
        J_x_low = zeros(n_x);
        J_p_low = zeros(n_x,n_p);

        % Parameters
        M = [25.8 0 0; ...
             0 33.8 1.0115; ...
             0 1.0115 2.16];
        N = [2 0 0; ...
             0 7 0.1; ...
             0 0.1 0.5];
        k1 = 25;
        k2 = 2.5;
        k3 = 5;

        inv_det_M23 = 1/det(M(2:3,2:3));

        % Constant elements
        J_x_low(3,6) = 1;
        J_x_low(4,4) = -k1/M(1,1);
        J_x_low(5,3) = k2*M(2,3)*inv_det_M23;
        J_x_low(5,5) = (N(3,2)*M(2,3) - N(2,2)*M(3,3))*inv_det_M23;
        J_x_low(5,6) = (-N(2,3)*M(3,3) + (k3+N(3,3))*M(2,3))*inv_det_M23;
        J_x_low(6,3) = -k2*M(2,2)*inv_det_M23;
        J_x_low(6,5) = (N(2,2)*M(3,2) - N(3,2)*M(2,2))*inv_det_M23;
        J_x_low(6,6) = (N(2,3)*M(3,2) - (k3+N(3,3))*M(2,2))*inv_det_M23;

        J_x_up = J_x_low;

        % Bounded elements (with all cos/sin between -1 and 1)
        J_x_low(1,3) = min([x_low(4),x_up(4),-x_low(4),-x_up(4)]) + min([x_low(5),x_up(5),-x_low(5),-x_up(5)]);
        J_x_low(2,3) = min([x_low(4),x_up(4),-x_low(4),-x_up(4)]) + min([x_low(5),x_up(5),-x_low(5),-x_up(5)]);
        J_x_low(1:2,4:5) = -1;
        
        J_x_up(1,3) = max([x_low(4),x_up(4),-x_low(4),-x_up(4)]) + max([x_low(5),x_up(5),-x_low(5),-x_up(5)]);
        J_x_up(2,3) = max([x_low(4),x_up(4),-x_low(4),-x_up(4)]) + max([x_low(5),x_up(5),-x_low(5),-x_up(5)]);
        J_x_up(1:2,4:5) = 1;
        
        % Input Jacobian (all elements are constant)
        J_p_low(4,1) = k1/M(1,1);
        J_p_low(5:6,2) = M(2:3,2:3)\[0;k2];
        J_p_up = J_p_low;
        
    case 4
        %% Traffic diverge 3 links (discrete-time)
        v = 0.5;            % free-flow speed, in links/period
        w = 1/6;            % congestion-wave speed, in links/period

        J_x_low = [1-v 0 0;0 1-(v+w) -w;0 -w 1-(v+w)];
        J_x_up = [1 2*w 2*w;v/2 1 0;v/2 0 1];
        J_p_low = [1;0;0];     % input: d1
        J_p_up = J_p_low;
        
    case 5
        %% Temperature variations in circular building (discrete-time)
        % Parameters
        a = 0.45;
        b = 0.045;
        c = 0.09;
        Th = 50;    % Heater temperature
        
        n_x = length(x_low);
        J_x_low = (1-2*a-b)*eye(n_x) -c*diag(p_up) + a*circshift(eye(n_x),1) + a*circshift(eye(n_x),-1);
        J_x_up = (1-2*a-b)*eye(n_x) -c*diag(p_low) + a*circshift(eye(n_x),1) + a*circshift(eye(n_x),-1);
        
        J_p_low = c*(Th-diag(x_up));
        J_p_up = c*(Th-diag(x_low));

    case 9
        %% Traffic diverge nx-link (continuous-time ; needs n_x >= 5)
        % This system corresponds to a traffic flow like this:
        %         ---(x2)---(x4)--- ... ---(x(end-1))
        %        /
        % (x1)---
        %        \
        %         ---(x3)---(x5)--- ... ---(x(end))

        % Parameters
        v = 0.5;            % free-flow speed, in links/period
        w = 1/6;            % congestion-wave speed, in links/period
        b = 3/4;            % fraction of vehicule staying on the network after each link
        T = 30;             % time step for the continuous-time model, in seconds

        % State Jacobian
        J_x_low = zeros(n_x);
        J_x_up = zeros(n_x);
        
        J_x_low(1:3,1:5) = 1/T*[-v 0 0 0 0;0 -(v+w) -w 0 0;0 -w -(v+w) 0 0];
        J_x_up(1:3,1:5) = 1/T*[0 2*w 2*w 0 0;v/2 0 0 w/b 0;v/2 0 0 0 w/b];
        for i=4:(n_x-2)
            J_x_low(i,i) = -(w+v)/T;
            J_x_up(i,i-2) = b*v/T;
            J_x_up(i,i+2) = w/b/T;
        end
        J_x_low(end, end) = -(w+v)/T;
        J_x_low(end-1, end-1) = -(w+v)/T;
        J_x_up(end, end-2) = b*v/T;
        J_x_up(end-1, end-3) = b*v/T;

        % Input Jacobian
        J_p_low = zeros(n_x);
        J_p_low(1,1) = 1;
        J_p_up = J_p_low;
        
    case 10
        %% Nonlinear Longitudinal Model of a Quadrotor (continuous-time)
        J_x_low = zeros(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up = zeros(n_x,n_p);
        
        % Parameters
        K = 0.89/1.4;
        d0 = 70;
        d1 = 17;

        % Constant terms
        J_x_low(1,3) = 1;
        J_x_low(2,4) = 1;
        J_x_low(5,6) = 1;
        J_x_low(6,5:6) = [-d0,-d1];
        J_x_up = J_x_low;
        
        % Time varying terms to be bounded
        % (assuming that state constraints x(5) in [-pi/12,pi/12])
        x5min = -pi/12;
        x5max = pi/12;
        if u(1)>=0
            J_x_low(3,5) = u(1)*K*min(cos(x5min),cos(x5max));
            J_x_up(3,5) = u(1)*K;
            J_x_low(4,5) = -u(1)*K*sin(x5max);
            J_x_up(4,5) = -u(1)*K*sin(x5min);
        else
            J_x_low(3,5) = u(1)*K;
            J_x_up(3,5) = u(1)*K*min(cos(x5min),cos(x5max));
            J_x_low(4,5) = -u(1)*K*sin(x5min);
            J_x_up(4,5) = -u(1)*K*sin(x5max);
        end

    case 11
        %% Pursuer-evader game with 2 Dubin's vehicles (continuous-time)
        time_step = t_final-t_init; 
        J_x_low = zeros(n_x);
        J_x_up = zeros(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up = zeros(n_x,n_p);
         
        % Constant terms:
        J_p_low(6,2) = 1;
        J_p_up(6,2) = 1;

        % Local bounds on x(3) to bound the Jacobians
        % Possible orientation values over the next sampling period
        theta_low = x_low(3) + min(0,time_step*u(2));
        theta_up = x_up(3) + max(0,time_step*u(2));
        % Min/Max for cosine
        % (conditions only valid if x(3) initial set is in [-pi,pi])
        min_cos = -1;
        if ~((theta_low <= -pi && theta_up >= -pi) || (theta_low <= pi && theta_up >= pi))
            min_cos = min(cos(theta_low),cos(theta_up));
        end
        max_cos = 1;
        if ~((theta_low <= -2*pi && theta_up >= -2*pi) || (theta_low <= 0 && theta_up >= 0) || (theta_low <= 2*pi && theta_up >= 2*pi))
            max_cos = max(cos(theta_low),cos(theta_up));
        end
        % Min/Max for sine
        % (conditions only valid if x(3) initial set is in [-pi,pi])
        min_sin = -1;
        if ~((theta_low <= -5*pi/2 && theta_up >= -5*pi/2) || (theta_low <= -pi/2 && theta_up >= -pi/2) || (theta_low <= 3*pi/2 && theta_up >= 3*pi/2))
            min_sin = min(sin(theta_low),sin(theta_up));
        end
        max_sin = 1;
        if ~((theta_low <= -3*pi/2 && theta_up >= -3*pi/2) || (theta_low <= pi/2 && theta_up >= pi/2) || (theta_low <= 5*pi/2 && theta_up >= 5*pi/2))
            max_sin = max(sin(theta_low),sin(theta_up));
        end
        % Bounds of the partial derivatives
        if u(1) > 0
            J_x_low(1,3) = -u(1)*max_sin;
            J_x_up(1,3) = -u(1)*min_sin;
            J_x_low(2,3) = u(1)*min_cos;
            J_x_up(2,3) = u(1)*max_cos;
        else
            J_x_low(1,3) = -u(1)*min_sin;
            J_x_up(1,3) = -u(1)*max_sin;
            J_x_low(2,3) = u(1)*max_cos;
            J_x_up(2,3) = u(1)*min_cos; 
        end
        
        % Local bounds on x(6) to bound the Jacobians
        % Possible orientation values over the next sampling period
        theta_low = x_low(6) + min(0,time_step*p_low(2));
        theta_up = x_up(6) + max(0,time_step*p_up(2));
        % Min/Max for cosine
        % (conditions only valid if x(6) initial set is in [-pi,pi])
        min_cos = -1;
        if ~((theta_low <= -pi && theta_up >= -pi) || (theta_low <= pi && theta_up >= pi))
            min_cos = min(cos(theta_low),cos(theta_up));
        end
        max_cos = 1;
        if ~((theta_low <= -2*pi && theta_up >= -2*pi) || (theta_low <= 0 && theta_up >= 0) || (theta_low <= 2*pi && theta_up >= 2*pi))
            max_cos = max(cos(theta_low),cos(theta_up));
        end
        % Min/Max for sine
        % (conditions only valid if x(6) initial set is in [-pi,pi])
        min_sin = -1;
        if ~((theta_low <= -5*pi/2 && theta_up >= -5*pi/2) || (theta_low <= -pi/2 && theta_up >= -pi/2) || (theta_low <= 3*pi/2 && theta_up >= 3*pi/2))
            min_sin = min(sin(theta_low),sin(theta_up));
        end
        max_sin = 1;
        if ~((theta_low <= -3*pi/2 && theta_up >= -3*pi/2) || (theta_low <= pi/2 && theta_up >= pi/2) || (theta_low <= 5*pi/2 && theta_up >= 5*pi/2))
            max_sin = max(sin(theta_low),sin(theta_up));
        end
        % Bounds of the partial derivatives
        J_x_low(4,6) = min([-p_low(1)*min_sin,-p_low(1)*max_sin,-p_up(1)*min_sin,-p_up(1)*max_sin]);
        J_x_up(4,6) = max([-p_low(1)*min_sin,-p_low(1)*max_sin,-p_up(1)*min_sin,-p_up(1)*max_sin]);
        J_x_low(5,6) = min([p_low(1)*min_cos,p_low(1)*max_cos,p_up(1)*min_cos,p_up(1)*max_cos]);
        J_x_up(5,6) = max([p_low(1)*min_cos,p_low(1)*max_cos,p_up(1)*min_cos,p_up(1)*max_cos]);
        J_p_low(4,1) = min_cos; 
        J_p_up(4,1) = max_cos; 
        J_p_low(5,1) = min_sin; 
        J_p_up(5,1) = max_sin; 
        
    case 12
        %% Simple oscillator model (continuous-time)
        J_x_low = [0 1;-1 0];
        J_x_up = J_x_low;
        J_p_low = zeros(n_x,n_p);
        J_p_up = J_p_low;

end


