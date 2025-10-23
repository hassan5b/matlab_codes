%% Product a*B of scalar interval a with interval matrix B
% Input. [a_LB, a_UB]: scalar interval a
% Input. [B_LB, B_UB]: interval matrices B
% Outputs. [prod_LB, prod_UB]: interval matrix for the product a with each element of B

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [prod_LB,prod_UB] = Interval_matrix_product_with_scalar(a_LB,a_UB,B_LB,B_UB)
prod_LB = min(min(min(a_LB*B_LB,a_LB*B_UB),a_UB*B_LB),a_UB*B_UB);
prod_UB = max(max(max(a_LB*B_LB,a_LB*B_UB),a_UB*B_LB),a_UB*B_UB);
end
