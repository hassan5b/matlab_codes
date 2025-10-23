%% Reachable set over-approximation using quasi-Monte Carlo simulations for continuous or discrete-time systems
% This method is similar to a Monte Carlo simulation but uses a
% deterministic scheme for selecting sample points, allowing the method
% to guarantee an over-approximation.
% This method can also be seen as a generalization of method 2 in
% OA_methods/OA_2_CT_Contraction_growth_Bound.m to a grid sampling of the
% intervals of initial states and inputs, instead of a single central
% sample point as in method 2.

% Computation proceeds as:
%   1. Given a number of sample points to sample along each dimension of
%      the initial state interval and the input interval, set up a uniform
%      grid over the intervals, and take every pair of initial state
%      and input as a sample point.
%   2. For each sample point, simulate the system trajectory to generate
%      a successor. 
%   3. Using the error bound function, compute the maximum distance 
%      by which two samples, given the spacing in the grids over the 
%      initial state and input intervals, can differ in their successors.
%   4. Taking the interval hull of the successors, add this maximum error
%      to the upper bound, and subtract it from the lower bound.
%      This resulting interval is the reachable set over-approximation.

% This method guarantees an over-approximation, unlike the Monte Carlo
% method (OA_methods/OA_7_DT_CT_Monte_Carlo.m). However, the Monte Carlo
% method's time performance can scale better in the sizes of the input
% and state dimensions compared to quasi-Monte Carlo.

% Source paper:
% P.-J. Meyer, A. Devonport and M. Arcak, "Sampling-Based Methods", 
% In: "Interval Reachability Analysis: Bounding Trajectories of Uncertain 
% Systems with Boxes for Control and Verification". Springer, 2021
% DOI: 10.1007/978-3-030-65110-7_7

% List of inputs
%   time_vector: time information for the reachability analysis
%       1D vector for discrete-time system: time_vector = t_init
%       2D vector for continuous-time system: time_vector=[t_init, t_final]
%       t_init: initial time
%       t_final: if present, time at which reachable set is approximated.
%           For discrete-time systems, t_final will be one step from t_init
%   [x_low, x_up]: interval of initial states (at time t_init).
%   [p_low, p_up]: interval of allowed input values.
%   QMC_samples_per_state_dimension: vector of same length as state variable. 
%       The i'th value of this vector is the number of 
%       points to sample along dimension i of the initial state interval.
%   QMC_samples_per_input_dimension: vector of same length as input variable.
%       The j'th value of this vector is the number of points
%       to sample along dimension j of the interval over input values.
%   error_bound_handle: function handle of the error bound.
%       For discrete-time systems:
%           inputs: t_init, infinity norm error of state, infinity norm
%               error of input.
%           outputs: non-negative real error bound.
%       For continuous-time systems:
%           inputs: [t_init, t_final], infinity norm error of state,
%               infinity norm error of input.
%           outputs: non-negative real error bound.

% List of outputs
%   [succ_low, succ_up]: interval over-approximation of the reachable set
%       at time t_final for continuous-time systems,
%       or after one time step for discrete-time systems

% Authors:  
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [succ_low, succ_up] = OA_6_DT_CT_Quasi_Monte_Carlo(time_vector,...
    x_low, x_up, p_low, p_up,...
    QMC_samples_per_state_dimension, QMC_samples_per_input_dimension,...
    error_bound_handle)

% State and input dimensions
n_x = length(x_low);
n_p = length(p_low);

% Extract ODE solver to use.
run('Solver_parameters.m');
ode_solver = parameters.ode_solver;
ode_options = parameters.ode_options;

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

% Verify number of points per state dimension.
assert(length(QMC_samples_per_state_dimension) == length(x_low),...
    'size of QMC_samples_per_state_dimension must match size of state x_low');
assert(all(QMC_samples_per_state_dimension >= 1),...
    'The number of samples per state dimension must be a positive integer');
assert(all(mod(QMC_samples_per_state_dimension, 1) == 0),...
    'The number of samples per state dimension must be a positive integer');

% Verify number of points per input dimension.
assert(length(QMC_samples_per_input_dimension) == length(p_low),...
    'size of QMC_samples_per_input_dimension must match size of input p_low');
assert(all(QMC_samples_per_input_dimension >= 1),...
    'The number of samples per input dimension must be a positive integer');
assert(all(mod(QMC_samples_per_input_dimension, 1) == 0),...
    'The number of samples per input dimension must be a positive integer');

%% Initialization
state_range = x_up - x_low;
input_range = p_up - p_low;

% Grid spacing on initial state interval and input interval.
state_dimension_spacing = state_range ./ QMC_samples_per_state_dimension;
input_dimension_spacing = input_range ./ QMC_samples_per_input_dimension;

sample_size = prod(QMC_samples_per_state_dimension) * prod(QMC_samples_per_input_dimension);
successor_samples = zeros(sample_size, n_x);

combined_num_points = vertcat(QMC_samples_per_state_dimension, QMC_samples_per_input_dimension);

%% Grid-based sampling and evaluation

% Compute sample coordinates.
sample_coords = cell(n_x + n_p, 1);
[sample_coords{:}] = ind2sub(combined_num_points, 1:sample_size);
sample_coords = cell2mat(sample_coords) - 1;

for sample_index = 1:sample_size
    % Create sample point from sample_coords.
    initial_state_sample = x_low + state_dimension_spacing .* sample_coords(1:n_x, sample_index) + state_dimension_spacing/2;
    input_sample = p_low + input_dimension_spacing .* sample_coords(n_x+1:n_x+n_p, sample_index) + input_dimension_spacing/2;
    
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

%% Final over-approximation of the reachable set

% Compute error bounds required to obtain an over-approximation
error = error_bound_handle(time_vector, norm(state_dimension_spacing,Inf)/2, norm(input_dimension_spacing,Inf)/2);

% Reachable set over-approximation:
%   interval hull of the sample point evaluations
%   expanded by the error bound on both sides
succ_low = min(successor_samples)' - error;
succ_up = max(successor_samples)' + error;
