%% Power operator of an interval matrix A^n
% Input. [A_LB,A_UB]: interval matrix A
% Input. n: power integer
% Outputs. [power_LB,power_UB]: interval matrix for A^n

% Authors:  
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

function [power_LB,power_UB] = Interval_matrix_power(A_LB,A_UB,n)

import interval_arithmetics.*;

if  n == 1  % A^1 = A
    power_LB = A_LB;
    power_UB = A_UB;
else        % A^n = A^(n-1)*A
    [temp_LB,temp_UB] = Interval_matrix_power(A_LB,A_UB,n-1);
    [power_LB,power_UB] = Interval_matrix_product(temp_LB,temp_UB,A_LB,A_UB);
    
    % A*A^(n-1) and A^(n-1)*A yield different results for n>2,
    % but both are guaranteed to over-approximate the real value of A^n
    % To use A^n = A*A^(n-1), replace the last line by the following:
%     [power_LB,power_UB] = Interval_matrix_product(A_LB,A_UB,temp_LB,temp_UB);
end
end