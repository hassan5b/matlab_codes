%% Norm of width of interval matrix
% Input:
%   [A_low, A_up]: Interval matrix.
% Output:
%   width: infinity norm of A_up - A_low.

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function width = Interval_matrix_width_norm(A_low, A_up)
    width = norm(A_up - A_low, Inf);
end