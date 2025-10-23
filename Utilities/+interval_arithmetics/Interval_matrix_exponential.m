%% Over-approximation of exponential of interval matrix e^[A_low, A_up]
% Exponential of interval matrix using scaling and squaring,
% and the Horner formulation of the truncated Taylor series as the
% submethod for computing the scaled down exponential.

% Source paper:
% A. Goldsztejn and A. Neumaier, "On the exponentiation of interval 
% matrices". Reliable Computing, v. 20, pp. 53-72, 2014.

% Source paper 2:
% C. Moler and C. Van Loan, "Nineteen dubious ways to compute the 
% exponential of a matrix". SIAM review, v. 20, n. 4, pp. 801-836, 1978.

% List of inputs
%   [A_low, A_up]: interval matrix A.
%   taylor_order: truncated Taylor series of exponential order.
%   scaling_factor: Scaling factor in scaling and squaring exponentials method.

% List of outputs
%   [X_low, X_up]: interval matrix of exponential of [A_low, A_up]

% Authors:  
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [X_low, X_up] = Interval_matrix_exponential(A_low, A_up, taylor_order, scaling_factor)
    import interval_arithmetics.Interval_matrix_norm;
    import interval_arithmetics.Interval_matrix_square;

    %% Scaling and squaring process
    
    % Minimum scaling factor L such that ||A||/(2^L) <= 1
    %   for the results to be better than the classical method
    scaling_factor_min = ceil(log2(Interval_matrix_norm(A_low, A_up)));
    if (nargin < 4) || isnan(scaling_factor) || (scaling_factor < scaling_factor_min)
        scaling_factor = scaling_factor_min;
    end
    
    % Minimum taylor order to obtain a sound over-approximation:
    %   need (taylor_order + 2) * 2^(scaling_factor) > ||[A_low, A_up]||
    taylor_order = max(taylor_order,ceil(Interval_matrix_norm(A_low, A_up)/2^(scaling_factor)-2));
    
    % Horner formulation of the scaled-down exponential
    [X_low, X_up] = Horner_exponential(A_low / (2^scaling_factor), A_up / (2^scaling_factor), taylor_order);

    % Scale back up the obtained interval: X^(2^scaling_factor).
    for i = 1:scaling_factor
        [X_low, X_up] = Interval_matrix_square(X_low, X_up);
    end
end

%% Horner formulation of Taylor series of exponential.
% Achieves a tighther interval than the direct computation of
% the Taylor series.
% Inputs:
%   [A_low, A_up]: interval matrices A.
%   taylor_order: truncation order of the Taylor series.
% Outputs:
%   [X_low, X_up]: interval matrix of exponential of [A_low, A_up]
function [X_low, X_up] = Horner_exponential(A_low, A_up, taylor_order)
    import interval_arithmetics.Interval_matrix_product;
    import interval_arithmetics.Interval_matrix_sum;
    import interval_arithmetics.Interval_matrix_norm;
    
    assert(taylor_order + 2 > Interval_matrix_norm(A_low, A_up),...
        'taylor order + 2 > ||[A_low, A_up]|| must be satisfied for Horner exponential to return an overapproximation');
    
    n_x = size(A_low, 1);

    [X_low, X_up] = Interval_matrix_sum(eye(n_x), eye(n_x),...
        A_low / taylor_order, A_up / taylor_order);
    for i = (taylor_order-1):-1:1
        [X_low, X_up] = Interval_matrix_product(A_low / i, A_up / i, X_low, X_up);
        [X_low, X_up] = Interval_matrix_sum(eye(n_x), eye(n_x), X_low, X_up);
    end
    
    % Add in over-approximation of residual terms.
    [E_low, E_up] = Truncation_remainder(A_low, A_up, taylor_order);
    [X_low, X_up] = Interval_matrix_sum(X_low, X_up, E_low, E_up);
end

%% Old implementation of the interval matrix exponential
% More straightforward definition, but also more conservative
% Inputs:
%   [A_low, A_up]: interval matrices A.
%   taylor_order: truncation order of the Taylor series.
% Outputs:
%   [X_low, X_up]: interval matrix of exponential of [A_low, A_up]
function [X_low, X_up] = Taylor_exponential(A_low, A_up, taylor_order)
    import interval_arithmetics.Interval_matrix_product;
    import interval_arithmetics.Interval_matrix_sum;
    import interval_arithmetics.Interval_matrix_norm;
    
    assert(taylor_order + 2 > Interval_matrix_norm(A_low, A_up),...
        'taylor order + 2 > ||[A_low, A_up]|| must be satisfied for Taylor exponential to return an overapproximation');
    
    n_x = size(A_low, 1);
    X_low = eye(n_x);
    X_up = eye(n_x);
    
    new_term_low = eye(n_x);
    new_term_up = eye(n_x);
    
    for i = 1:taylor_order
        % Compute A^i/i! based on term of i-1
        [new_term_low, new_term_up] = Interval_matrix_product(...
            new_term_low, new_term_up, A_low / i, A_up / i);
        [X_low, X_up] = Interval_matrix_sum(X_low, X_up, new_term_low, new_term_up);
    end
    
    % Add in over-approximation of residual terms.
    [E_low, E_up] = Truncation_remainder(A_low, A_up, taylor_order);
    [X_low, X_up] = Interval_matrix_sum(X_low, X_up, E_low, E_up);
end

%% Over-approximation of the remainder of the truncated Taylor series
% Inputs:
%   [A_low, A_up]: interval matrices A.
%   taylor_order: truncation order of the Taylor series.
% Outputs:
%   [E_low, E_up]: bounds of the remainder of the truncated Taylor series
function [E_low,E_up] = Truncation_remainder(A_low, A_up, taylor_order)
    import interval_arithmetics.Interval_matrix_norm;

    [n_x,~] = size(A_low);
    A_norm = Interval_matrix_norm(A_low,A_up);
    E_up = (A_norm^(taylor_order + 1)) / (factorial(taylor_order + 1) * (1 - A_norm / (taylor_order + 2)));
    E_up = ones(n_x)*E_up;
    E_low = -E_up;
end
