%% Sum of two interval matrices A+B
% Inputs. [A_LB, A_UB], [B_LB, B_UB]: two interval matrices A and B
% Outputs. [sum_LB, sum_UB]: interval matrix for the sum A+B

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [sum_LB, sum_UB] = Interval_matrix_sum(A_LB, A_UB, B_LB, B_UB)
sum_LB = A_LB + B_LB;
sum_UB = A_UB + B_UB;
end