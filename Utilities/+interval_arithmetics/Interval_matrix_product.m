%% Product of two interval matrices A*B
% Inputs. [A_LB,A_UB], [B_LB,B_UB]: two interval matrices A and B
% Outputs. [prod_LB,prod_UB]: interval matrix for the product A*B

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [prod_LB,prod_UB] = Interval_matrix_product(A_LB,A_UB,B_LB,B_UB)
prod_LB = zeros(size(A_LB,1),size(B_LB,2));
prod_UB = zeros(size(A_LB,1),size(B_LB,2));
for i = 1:size(A_LB,1)  % lines of A
    for j = 1:size(B_LB,2)  % columns of B
        prod_LB(i,j) = sum(min([A_LB(i,:)'.*B_LB(:,j) A_UB(i,:)'.*B_LB(:,j) A_LB(i,:)'.*B_UB(:,j) A_UB(i,:)'.*B_UB(:,j)],[],2));
        prod_UB(i,j) = sum(max([A_LB(i,:)'.*B_LB(:,j) A_UB(i,:)'.*B_LB(:,j) A_LB(i,:)'.*B_UB(:,j) A_UB(i,:)'.*B_UB(:,j)],[],2));
    end
end
end