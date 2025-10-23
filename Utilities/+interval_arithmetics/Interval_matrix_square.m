%% Square of an interval matrix
% Computes the optimal square of an interval matrix, removing interval
% dependency.

% Source Paper:
% O. Koshaleva, V. Kreinovich, G. Mayer, and H. T. Nguyen. "Computing the
% Cube of an Interval Matrix is NP-hard". Departmental Technical Reports
% (CS), 2005.

% Input. [A_low, A_up]: nxn interval matrix A
% Outputs. [square_low, square_up]: square of nxn interval matrix A

% Authors:  
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [square_low, square_up] = Interval_matrix_square(A_low, A_up)
    import interval_arithmetics.*;
    [n1, n2] = size(A_low);
    assert(n1 == n2, 'Matrix must be square');
    assert(n1 == size(A_up, 1), 'A_low must be same size as A_up');
    assert(n2 == size(A_up, 2), 'A_low must be same size as A_up');
        
    n = n1;
    
    square_low = zeros(n);
    square_up  = zeros(n);
    
    for i = 1:n
        for j = 1:n
            if i ~= j
                for k = 1:n
                    if k == i || k == j
                        continue;
                    end
                    [c_low, c_up] = Interval_matrix_product(A_low(i, k), A_up(i, k), A_low(k, j), A_up(k, j));
                    [square_low(i, j), square_up(i, j)] = Interval_matrix_sum(square_low(i, j), square_up(i, j), c_low, c_up);
                end
                [c_low, c_up] = Interval_matrix_sum(A_low(i, i), A_up(i, i), A_low(j, j), A_up(j, j));
                [c_low, c_up] = Interval_matrix_product(A_low(i, j), A_up(i, j), c_low, c_up);
                [square_low(i, j), square_up(i, j)] = Interval_matrix_sum(square_low(i, j), square_up(i, j), c_low, c_up);
            else
                assert(i == j);
                for k = 1:n
                    if k == i
                        [c_low, c_up] = Interval_scalar_square(A_low(i, i), A_up(i, i));
                    else
                        [c_low, c_up] = Interval_matrix_product(A_low(i, k), A_up(i, k), A_low(k, i), A_up(k, i));
                    end
                    [square_low(i, i), square_up(i, i)] = Interval_matrix_sum(square_low(i, i), square_up(i, i), c_low, c_up);
                end
            end                
        end
    end
end

% Optimal scalar interval square.
function [c_low, c_up] = Interval_scalar_square(a_low, a_up)
    assert(a_low <= a_up);
    if 0 <= a_low
        c_low = a_low^2;
        c_up = a_up^2;
    elseif a_up <= 0
        c_low = a_up^2;
        c_up = a_low^2;
    else
        c_low = 0;
        c_up = max(a_low^2, a_up^2);
    end
end