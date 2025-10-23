%% Infinity norm definition for an interval matrix
% Inputs. [A_LB,A_UB]: interval matrix A
% Output. inf_norm: infinity norm of interval matrix A

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function inf_norm = Interval_matrix_norm(A_LB,A_UB)
inf_norm = norm(max(abs(A_LB),abs(A_UB)),Inf);
end