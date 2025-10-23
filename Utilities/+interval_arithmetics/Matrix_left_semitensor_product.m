%% Left Semi-Tensor Product of Matrices using Kronecker product
% Left Semi-Tensor Product(A,B) is to calculate the (left) semi-tensor
% product of two matrices A and B. The number of columns of A must be
% either a divisor or a multiple of the number of rows of B.

% Source paper:
% D. Cheng, H. Qi and Y. Zhao, "An introduction to semi-tensor product of 
% matrices and its applications". World Scientific, 2012

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function C = Matrix_left_semitensor_product(A,B)
log = logger.get_logger();

if ~(isa(A,'sym') || isa(A,'double'))
    A = double(A);
end
if ~(isa(B,'sym') || isa(B,'double'))
    B = double(B);
end

assert((ndims(A)<=2) && (ndims(B)<=2),'Input arguments must be 2-D.')

n = size(A,2); % get the number of columns
p = size(B,1); % get the number of rows
if n == p
    C = A*B;
elseif mod(n,p) == 0
    C = A*kron(B,eye(n/p));
elseif mod(p,n) == 0
    C = kron(A,eye(p/n))*B;
else
    log.error('The input arguments do not meet the multiple dimension matching condition.')
end
end
