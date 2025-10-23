%% User-provided signs on the Jacobian matrices for continuous-time or discrete-time systems
% Used in:
% - OA_methods/TIRA.m
%       to check the requirements for over-approximation methods 1,3,5 in
%       OA_methods/OA_1_CT_Monotonicity.m
%       OA_methods/OA_3_CT_Mixed_Monotonicity.m
%       OA_methods/OA_5_DT_Mixed_Monotonicity.m

% The user can either provider global signs, or a function of the input
% arguments (t_init,t_final,x_low,x_up,p_low,p_up) returning local signs

% Jacobian definitions:
%   to states:  J_x(t) = d(System_description(t,x,p))/dx
%   to inputs:  J_p(t) = d(System_description(t,x,p))/dp

% List of inputs
%   t_init: initial time
%   t_final: time at which the reachable set is approximated (for continuous-time system only)
%       for discrete-time system, a dummy value can be provided
%   [x_low,x_up]: interval of initial states (at time t_init)
%   [p_low,p_up]: interval of allowed input values

% List of outputs
%   J_x_signs: sign matrix of the Jacobian with respect to the state
%   J_p_signs: sign matrix of the Jacobian with respect to the input

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function [J_x_signs,J_p_signs] = UP_Jacobian_Signs(t_init,t_final,x_low,x_up,p_low,p_up)
n_x = length(x_low);
n_p  = length(p_low);

%% Default values as NaN (not a number)
J_x_signs = NaN(n_x);
J_p_signs = NaN(n_x,n_p);

%% User-provided Jacobian signs
% Can be either global signs
% or local signs depending on inputs: t_init,t_final,x_low,x_up,p_low,p_up

% If System_description.m has no input variable 'p', uncomment below:
% J_p_signs = zeros(n_x,n_p);

global system_choice

switch system_choice
    case 4
        %% Traffic diverge 3 links (discrete-time)
        J_x_signs = [1 1 1;1 1 -1;1 -1 1];
        J_p_signs = [1;0;0];
        
    case 5
        %% Temperature variations in circular building (discrete-time)
        % These signs are only valid if heater_max<=0.6 in MAIN_CALL.m
        J_x_signs = eye(n_x)+circshift(eye(n_x),1)+circshift(eye(n_x),-1);
        J_p_signs = eye(n_x);
        
    case 8
        %% Protein interaction (continuous-time)
        J_x_signs = [-1 1;1 -1];
        J_p_signs = zeros(n_x,n_p);
        
    case 12
        %% Simple oscillator model (continuous-time)
        J_x_signs = [0 1;-1 0];
        J_p_signs = zeros(n_x,n_p);
        
end




