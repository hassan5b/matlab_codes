%% User-provided bounds on the second-order Jacobian matrices for a continuous-time system
% Used in:
% - Utilities/Sensitivity_bounds_choice.m
%       to compute bounds on the sensitivity matrices of a continuous-time
%       system, to be used in over-approximation method 4 in
%       OA_methods/OA_4_CT_Sampled_data_MM.m
%       with submethod (4-3)

% The user can either provider global bounds, or a function of the input
% arguments (t_init,t_final,x_low,x_up,p_low,p_up) returning local bounds

% Given the first-order Jacobian definitions:
%   J_x(t,x,p) = d(System_description(t,x,p))/dx
%   J_p(t,x,p) = d(System_description(t,x,p))/dp
% the second-order Jacobian are defined by:
%   J_xx(t,x,p) = dJ_x(t,x,p)/dx
%   J_xp(t,x,p) = dJ_x(t,x,p)/dp
%   J_pp(t,x,p) = dJ_p(t,x,p)/dp
%   J_px(t,x,p) = dJ_p(t,x,p)/dx

% List of inputs
%   t_init: initial time
%   t_final: time at which the reachable set is approximated
%   [x_low,x_up]: interval of initial states (at time t_init)
%   [p_low,p_up]: interval of allowed input values

% List of outputs
%   bounds of the 4 second-order Jacobians described as 2D matrices
%       [J_xx_low,J_xx_up]: of dimensions (n_x,n_x^2) 
%       [J_xp_low,J_xp_up]: of dimensions (n_x,n_x*n_p) 
%       [J_pp_low,J_pp_up]: of dimensions (n_x,n_p^2) 
%       [J_px_low,J_px_up]: of dimensions (n_x,n_x*n_p) 

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [J_xx_low,J_xx_up,J_xp_low,J_xp_up,J_pp_low,J_pp_up,J_px_low,J_px_up] = UP_Jacobian_2ndOrder_Bounds(t_init,t_final,x_low,x_up,p_low,p_up)
n_x = length(x_low);
n_p  = length(p_low);

%% Default values as NaN (not a number)
J_xx_low = NaN(n_x,n_x,n_x);
J_xp_low = NaN(n_x,n_x,n_p);
J_pp_low = NaN(n_x,n_p,n_p);
J_px_low = NaN(n_x,n_p,n_x);
J_xx_up = NaN(n_x,n_x,n_x);
J_xp_up = NaN(n_x,n_x,n_p);
J_pp_up = NaN(n_x,n_p,n_p);
J_px_up = NaN(n_x,n_p,n_x);

%% User-provided second-order Jacobian bounds as 3D matrices, with:
%   size(J_xx_low) = [n_x,n_x,n_x] and J_xx_low(i,j,k) = dJ_x(i,j)/d_x(k)
%   size(J_xp_low) = [n_x,n_x,n_p] and J_xp_low(i,j,k) = dJ_x(i,j)/d_p(k)
%   size(J_pp_low) = [n_x,n_p,n_p] and J_pp_low(i,j,k) = dJ_p(i,j)/d_p(k)
%   size(J_px_low) = [n_x,n_p,n_x] and J_px_low(i,j,k) = dJ_p(i,j)/d_x(k)
% These 3D matrices will be automatically converted to the appropriate 2D 
%   form at the end of this file.

% If System_description.m has no input variable 'p', uncomment below:
% J_xp_low = zeros(n_x,n_x,n_p);
% J_pp_low = zeros(n_x,n_p,n_p);
% J_px_low = zeros(n_x,n_p,n_x);
% J_xp_up = zeros(n_x,n_x,n_p);
% J_pp_up = zeros(n_x,n_p,n_p);
% J_px_up = zeros(n_x,n_p,n_x);

global system_choice
global u        % Some systems need an external control input

switch system_choice
    case 1
        %% Unicycle (continuous-time)
        J_xx_low = zeros(n_x,n_x,n_x);
        J_xp_low = zeros(n_x,n_x,n_p);
        J_pp_low = zeros(n_x,n_p,n_p);
        J_px_low = zeros(n_x,n_p,n_x);
        J_xx_up = zeros(n_x,n_x,n_x);
        J_xp_up = zeros(n_x,n_x,n_p);
        J_pp_up = zeros(n_x,n_p,n_p);
        J_px_up = zeros(n_x,n_p,n_x);
        
        time_step = t_final-t_init;

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
            J_xx_low(1,3,3) = -u(1)*max_cos;
            J_xx_up(1,3,3) = -u(1)*min_cos;
            J_xx_low(2,3,3) = -u(1)*max_sin;
            J_xx_up(2,3,3) = -u(1)*min_sin;
        else
            J_xx_low(1,3,3) = -u(1)*min_cos;
            J_xx_up(1,3,3) = -u(1)*max_cos;
            J_xx_low(2,3,3) = -u(1)*min_sin;
            J_xx_up(2,3,3) = -u(1)*max_sin;
        end
        
        J_xx_low(1,3,3) = -abs(u(1));
        J_xx_up(1,3,3) = abs(u(1));
        J_xx_low(2,3,3) = -abs(u(1));
        J_xx_up(2,3,3) = abs(u(1));
        
    case 2
        %% Traffic diverge 3 links (continuous-time)
        J_xx_low = zeros(n_x,n_x,n_x);
        J_xp_low = zeros(n_x,n_x,n_p);
        J_pp_low = zeros(n_x,n_p,n_p);
        J_px_low = zeros(n_x,n_p,n_x);
        J_xx_up = zeros(n_x,n_x,n_x);
        J_xp_up = zeros(n_x,n_x,n_p);
        J_pp_up = zeros(n_x,n_p,n_p);
        J_px_up = zeros(n_x,n_p,n_x);
        
    case 3
        %% Ship (continuous-time)
        J_xx_low = zeros(n_x,n_x,n_x);
        J_xp_low = zeros(n_x,n_x,n_p);
        J_pp_low = zeros(n_x,n_p,n_p);
        J_px_low = zeros(n_x,n_p,n_x);
        J_xx_up = zeros(n_x,n_x,n_x);
        J_xp_up = zeros(n_x,n_x,n_p);
        J_pp_up = zeros(n_x,n_p,n_p);
        J_px_up = zeros(n_x,n_p,n_x);

        % Ship (continuous-time): bounding all cos/sin between -1 and 1
        % J_xx
        J_xx_low(1:2,3,3) = min([x_low(4),-x_up(4)]) + min([x_low(5),-x_up(5)]);
        J_xx_low(1:2,3,4:5) = -1;
        J_xx_low(1:2,4:5,3) = -1;
        J_xx_up = -J_xx_low;
        % J_xp, J_pp, J_px: full zero matrices
          
    case 12
        %% Simple oscillator model (continuous-time)
        J_xx_low = zeros(n_x,n_x,n_x);
        J_xp_low = zeros(n_x,n_x,n_p);
        J_pp_low = zeros(n_x,n_p,n_p);
        J_px_low = zeros(n_x,n_p,n_x);
        J_xx_up = zeros(n_x,n_x,n_x);
        J_xp_up = zeros(n_x,n_x,n_p);
        J_pp_up = zeros(n_x,n_p,n_p);
        J_px_up = zeros(n_x,n_p,n_x);
               
end

%% Convert these 3D matrices to their corresponding 2D form:
%   size(J_xx_low) = [n_x,n_x^2]
%   size(J_xp_low) = [n_x,n_x*n_p]
%   size(J_pp_low) = [n_x,n_p^2]
%   size(J_px_low) = [n_x,n_x*n_p]
% and for example:
%   J_xx_low = [J_xx_low(:,1,:) , ... , J_xx_low(:,n_x,:)]
J_xx_low = reshape(permute(J_xx_low,[1,3,2]),[n_x,n_x^2]);
J_xx_up = reshape(permute(J_xx_up,[1,3,2]),[n_x,n_x^2]);
J_xp_low = reshape(permute(J_xp_low,[1,3,2]),[n_x,n_x*n_p]);
J_xp_up = reshape(permute(J_xp_up,[1,3,2]),[n_x,n_x*n_p]);
J_pp_low = reshape(permute(J_pp_low,[1,3,2]),[n_x,n_p^2]);
J_pp_up = reshape(permute(J_pp_up,[1,3,2]),[n_x,n_p^2]);
J_px_low = reshape(permute(J_px_low,[1,3,2]),[n_x,n_x*n_p]);
J_px_up = reshape(permute(J_px_up,[1,3,2]),[n_x,n_x*n_p]);
