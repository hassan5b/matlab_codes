%% User-provided description of the Jacobian functions
% Used in:
% - Utilities/Sensitivity_bounds_choice.m and
%   Utilities/Sensitivity_solver.m
%       when the solver of the sensitivity equations (for submethods (4-0) 
%       and (4-3) of over-approximation method 4 in 
%       OA_methods/OA_4_CT_Sampled_data_MM.m) is picked in 
%       Solver-parameters.m as the ODE solver on the sensitivity ODE
%       defined using explicit Jacobian function 
%       (sensitivity_equation_solver=2)

% Jacobian definitions:
%   to states:  J_x(t) = d(System_description(t,x,p))/dx
%   to inputs:  J_p(t) = d(System_description(t,x,p))/dp

% List of inputs
%   t_init: initial time
%   t: current time (for continuous-time system only)
%       for discrete-time system, a dummy value can be provided
%   x: current state
%   p: current input

% List of outputs
%   J_x: evaluation of the current Jacobian with respect to the state
%   J_p: evaluation of the current Jacobian with respect to the input

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function [J_x,J_p] = UP_Jacobian_Function(t_init,t,x,p)
n_x = length(x);
n_p  = length(p);

%% Default values as NaN (not a number)
J_x = NaN(n_x);
J_p = NaN(n_x,n_p);

%% User-provided Jacobian functions depending on inputs: t_init,t,x,p

% If System_description.m has no input variable 'p', uncomment below:
% J_p = zeros(n_x,n_p);



