function [P, L, U, sign] = splu(A)

% splu  Square PA=LU factorization *with row exchanges*.
%
% [P, L, U] = splu(A), for a square, invertible matrix A,
% uses Gaussian elimination to compute a permutation
% matrix P, a lower triangular matrix L and 
% an upper triangular matrix U so that P*A = L*U.
% P, L and U are the same size as A.
% sign = det(P); it is 1 or -1.
%
% See also slu, lu, rref, partic, nulbasis, colbasis.

[m, n] = size(A);
if m ~= n
   error('Matrix must be square.')
end
P = eye(n, n);
L = eye(n, n);
U = zeros(n, n);
tol = sqrt(eps);
sign = 1;

for k = 1:n
   if abs(A(k, k)) < tol
      for r = k:n
         if abs(A(r, k)) >= tol
            break
         end
         if r == n
            if nargout == 4
               sign = 0;
               return
            else
               disp('A is singular within tolerance.')
               error(['No pivot in column ' int2str(k) '.'])
            end
         end
      end
      A([r k], 1:n) = A([k r], 1:n);
      if k > 1, L([r k], 1:k-1) = L([k r], 1:k-1); end
      P([r k], 1:n) = P([k r], 1:n);
      sign = -sign;
   end
   for i = k+1:n
      L(i, k) = A(i, k) / A(k, k);
      for j = k+1:n
         A(i, j) = A(i, j) - L(i, k)*A(k, j);
      end
   end
   for j = k:n
      U(k, j) = A(k, j) * (abs(A(k, j)) >= tol);
   end
end

if nargout < 4
   roworder = P*(1:n)';
   disp('Pivots in rows:'), disp(roworder'); end
end
