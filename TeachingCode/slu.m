function [L, U] = slu(A)

% slu  LU factorization of a square matrix using *no row exchanges*.
%
% [L, U] = slu(A) uses Gaussian elimination to compute a unit 
% lower triangular L and an upper triangular U so that L*U = A.
% The algorithm will stop if a pivot entry is very small.
%
% See also slv, plu, lu.

[n, n] = size(A);

for k = 1:n
   if abs(A(k, k)) < sqrt(eps)
      disp(['Small pivot encountered in column ' int2str(k) '.'])
   end
   L(k, k) = 1;
   for i = k+1:n
      L(i,k) = A(i, k) / A(k, k);
      for j = k+1:n
         A(i, j) = A(i, j) - L(i, k)*A(k, j);
      end
   end
   for j = k:n
      U(k, j) = A(k, j);
   end
end
