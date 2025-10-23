%% Left Semi-Tensor Product of two interval matrices A*B

% Source paper:
% D. Cheng, H. Qi and Y. Zhao, "An introduction to semi-tensor product of 
% matrices and its applications". World Scientific, 2012

% Inputs. [A_LB,A_UB], [B_LB,B_UB]: two interval matrices A and B
%   The number of columns of A must be either a divisor or a multiple 
%   of the number of rows of B.

% Outputs. [prod_LB,prod_UB]: interval matrix for the product A*B
%   with the Left Semi-Tensor Product (using its Kronecker definition)

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [prod_LB,prod_UB] = Interval_left_semitensor_product(A_LB,A_UB,B_LB,B_UB)
log = logger.get_logger();
import interval_arithmetics.*;

% Check input size consistency
assert(isequal(size(A_LB),size(A_UB)),'The lower and upper bounds of the first interval matrix do not have the same dimensions')
assert(isequal(size(B_LB),size(B_UB)),'The lower and upper bounds of the second interval matrix do not have the same dimensions')
assert((ndims(A_LB)<=2) && (ndims(B_LB)<=2),'Input arguments must be 2-D matrices.')

n = size(A_LB,2); % get the number of columns
p = size(B_LB,1); % get the number of rows
if n == p
    [prod_LB,prod_UB] = Interval_matrix_product(A_LB,A_UB,B_LB,B_UB);
elseif mod(n,p) == 0
    z = n/p;
    [prod_LB,prod_UB] = Interval_matrix_product(A_LB,A_UB,kron(B_LB,eye(z)),kron(B_UB,eye(z)));
elseif mod(p,n) == 0
    z = p/n;
    [prod_LB,prod_UB] = Interval_matrix_product(kron(A_LB,eye(z)),kron(A_UB,eye(z)),B_LB,B_UB);
else
    log.error('The input arguments do not meet the multiple dimension matching condition.')
end
end