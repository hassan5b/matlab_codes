%% Over-approximation of the sensitivity bounds using interval arithmetics
% Requires user-provided Jacobian bounds in Input_files/UP_Jacobian_Bounds

% Source paper 1 (interval arithmetics for time-varying linear systems):
% M. Althoff, O. Stursberg and M. Buss, "Reachability analysis of linear 
% systems with uncertain parameters and inputs". 46th IEEE Conference on 
% Decision and Control, pp. 726-732, 2007. DOI: 10.1109/CDC.2007.4434084

% Source paper 2 (application of paper 1 to the sensitivity systems):
% P.-J. Meyer, S. Coogan and M. Arcak, "Sampled-data reachability analysis 
% using sensitivity and mixed-monotonicity". IEEE Control Systems Letters,
% v. 2(4), pp. 761-766, 2018. DOI: 10.1109/LCSYS.2018.2848280

% Source paper 3 (definition of basic operations for interval arithmetics):
% L. Jaulin, M. Kieffer, O. Didrit and E. Walter, "Applied interval 
% analysis: with examples in parameter and state estimation, robust control
% and robotics". Springer Science & Business Media, v. 1, 2001.

% List of inputs
%   t_init: initial time
%   t_final: time at which the reachable set is approximated
%   [x_low, x_up]: interval of initial states (at time t_init)
%   [p_low, p_up]: interval of allowed input values
%   [J_x_low, J_x_up]: bounds of the Jacobian with respect to the state
%   [J_p_low, J_p_up]: bounds of the Jacobian with respect to the input

% List of outputs
%   [S_x_low, S_x_up]: bounds of the sensitivity of successors (at time t_final) to variations of initial states (at time t_init)
%   [S_p_low, S_p_up]: bounds of the sensitivity of successors (at time t_final) to variations of inputs

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 2nd of December 2021

function [S_x_low,S_x_up,S_p_low,S_p_up] = Sensitivity_interval_arithmetics(t_init,t_final,x_low,x_up,p_low,p_up,J_x_low,J_x_up,J_p_low,J_p_up)

% State and input dimensions
n_x = length(x_low);
n_p = length(p_low);

import interval_arithmetics.*;
log = logger.get_logger();

% Check consistency of user-provided intervals
assert(all(J_x_low(:) <= J_x_up(:)),'Inconsistent bounds of the state Jacobian provided: need J_x_low <= J_x_up')
assert(all(J_p_low(:) <= J_p_up(:)),'Inconsistent bounds of the input Jacobian provided: need J_p_low <= J_p_up')

tic_IA = tic;

%% Computes the interval matrix over-approximating the sensitivity matrix with respect to initial conditions
if isequal(x_low, x_up)
    % No computation is needed if [x_low,x_up] contains a single value
    S_x_low = zeros(n_x);
    S_x_up  = zeros(n_x);
else
    % Use interval arithmetics to solve the set-valued system:
    % dS_x/dt = J_x*S_x, initialized with the identity matrix.
    % See source paper 2 for the problem definition and result statement
    % and source paper 1 for the technical details of the computations below.
    [S_x_low,S_x_up] = Affine_interval_system_reachability(t_init, t_final, eye(n_x), eye(n_x), J_x_low, J_x_up, zeros(n_x), zeros(n_x), 0);
end

%% Computes the interval matrix over-approximating the sensitivity matrix with respect to constant inputs
if isequal(p_low, p_up)
    % No computation is needed if [p_low,p_up] contains a single value
    S_p_low = zeros(n_x, n_p);
    S_p_up = zeros(n_x, n_p);
else
    % Use interval arithmetics to solve the set-valued system: 
    % dS_p/dt = J_x*S_p + J_p, initialized with the zero matrix.
    % See source paper 2 for the problem definition and result statement
    % and source paper 1 for the technical details of the computations below.
    [S_p_low, S_p_up] = Affine_interval_system_reachability(t_init, t_final, zeros(n_x, n_p), zeros(n_x, n_p), J_x_low, J_x_up, J_p_low, J_p_up, 0);
end

toc_IA = toc(tic_IA);
log.runtime('Computing sensitivity bounds took %f seconds\n', toc_IA);

