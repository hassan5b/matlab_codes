%% Reachable set approximation for any continuous or discrete-time system using Monte Carlo method
% Computation proceeds in three steps:
%   1. Given an accuracy parameter epsilon and confidence parameter delta,
%      compute the sample size N, according to equation (11) of
%      Devonport and Arcak (2020).
%   2. Select N initial states and N inputs uniformly at random from the
%      intervals [x_low, x_up] and [p_low, p_up] respectively. For each of
%      the N pairs, simulate the system trajectory with that initial state
%      and input, yielding N sample successor points that lie in the
%      reachable set.
%   3. Take the interval hull of these N successors to obtain an
%       approximation of the reachable set.

% The reachable set approximation produced by the Monte Carlo method is
% not guaranteed to be an over-approximation; however, it satisfies a
% probabilistic guarantee of accuracy and confidence. The details are in 
% (Devonport and Arcak (2020), Theorem 2), but here is a simple example.
% Suppose we have just computed a Monte Carlo estimate of the reachable
% set, and we draw a new (initial state, input) pair, that is, a pair that
% wasn't used to compute the estimate. What is the probability that this
% new sample yields a successor that is *not* in the Monte Carlo reachable
% set approximation we computed? With confidence (1-delta), this
% probability is no greater than epsilon, where the epsilon and delta are
% the parameters we used to compute the sample size.

% To ensure an over-approximation is generated, the quasi-Monte Carlo
% method can be applied: OA_methods/OA_6_DT_CT_Quasi_Monte_Carlo.m.
% However, the quasi-Monte Carlo method requires a number of samples that
% scales exponentially with the size of the state and input variables.

% This method has time complexity which is linear in the number of samples
% needed to satisfy the epsilon and delta parameters. The number of samples
% grows linearly with 1/epsilon, logarithmically with 1/delta, and linearly
% with nx (the size of the state variables). Therefore, delta can be 
% decreased (i.e. the confidence can be increased) with much less 
% performance penalty than can accuracy be increased by reducing epsilon.

% Source paper:
% Devonport, A. & Arcak, M.. (2020). Estimating Reachable Sets with
% Scenario Optimization. Proceedings of the 2nd Conference on Learning for
% Dynamics and Control, in PMLR 120:75-84 
 
% List of inputs
%   time_vector: time information for the reachability analysis
%       1D vector for discrete-time system: time_vector = t_init
%       2D vector for continuous-time system: time_vector=[t_init, t_final]
%       t_init: initial time
%       t_final: if present, time at which  reachable set is approximated.
%           For discrete-time systems, t_final will be one step from t_init
%   [x_low, x_up]: interval of initial states (at time t_init)
%   [p_low, p_up]: interval of allowed input values
%   epsilon: 1-epsilon is the accuracy
%   delta: 1-delta is the confidence in getting a (1-epsilon) accuracy 

% List of outputs
%   [succ_low, succ_up]: interval approximation of the reachable set
%       at time t_final for continuous-time systems,
%       or after one time step for discrete-time systems

% Authors:  
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [succ_low, succ_up] = OA_7_DT_CT_Monte_Carlo(time_vector, x_low, x_up, p_low, p_up, epsilon, delta)

% State and input dimensions
n_x = length(x_low);
n_p = length(p_low);

% Extract ODE solver to use.
run('Solver_parameters.m');
ode_solver = parameters.ode_solver;
ode_options = parameters.ode_options;
% Remove logger to be able to use the math function log()
clear log

%% Checking consistency of the input arguments

% Verify time vector size.
assert(length(time_vector) == 1 || length(time_vector) == 2,...
    'time_vector must be of length 1 for discrete-time, and of length 2 for continuous-time');

% Verify size and ordering of initial state interval.
assert(all(size(x_low) == size(x_up)),...
    'size of x_low must match size of x_up');
assert(all(x_up >= x_low),...
    'x_up must be greater than or equal to x_low, component-wise');

% Verify size and ordering of input interval.
assert(all(size(p_low) == size(p_up)),...
    'size of p_low must match size of p_up');
assert(all(p_up >= p_low),...
    'p_up must be greater than or equal to p_low, component-wise');

% Verify epsilon is scalar and in (0, 1).
assert(isscalar(epsilon), 'epsilon must be a scalar');
assert(epsilon > 0 && epsilon < 1, 'epsilon must be in (0, 1)');

% Verify delta is scalar and in (0, 1).
assert(isscalar(delta), 'delta must be a scalar');
assert(delta > 0 && delta < 1, 'delta must be in (0, 1)');

%% Initialization

% Determine number of samples needed to satisfy accuracy and confidence.
sample_size = ceil((exp(1)/(exp(1) - 1) * (log(1/delta) + 2 * n_x)) / epsilon);

% Array of successor samples
successor_samples = zeros(sample_size, n_x);

%% Random sampling and evaluation of the successors
for sample_index = 1:sample_size
    % Generate sample point.
    initial_state_sample = x_low + rand(n_x, 1) .* (x_up - x_low);
    input_sample = p_low + rand(n_p, 1) .* (p_up - p_low);
    
    % Evaluate the system successor for this sample point
    if length(time_vector) == 1
        % Discrete-time system
        successor_samples(sample_index, :) = System_description(time_vector, initial_state_sample, input_sample);
    else
        % Continuous-time system
        [~, successor_trajectory] = ode_solver(@(t, x) System_description(t, x, input_sample), time_vector, initial_state_sample, ode_options);
        successor_samples(sample_index, :) = successor_trajectory(end, :);
    end
end

%% Final approximation of the reachable set
%   as the interval hull of the computed successors
succ_low = min(successor_samples)';
succ_up = max(successor_samples)';
