function [P, L, U, pivcol, sign] = plu(A)

% plu  Rectangular PA=LU factorization *with row exchanges*.
%
% [P, L, U] = plu(A), for a rectangular matrix A, uses Gaussian elimination
% to compute a permutation matrix P, a lower triangular matrix L and 
% an upper trapezoidal matrix U so that PA = LU.
% U is the same size as A.  
% P and L are square, with as many rows as A.
% sign = det(P); it is 1 or -1.
%
% See also elim, slu, lu, rref, partic, nulbasis, colbasis.

[m, n] = size(A);
P = eye(m, m);
L = eye(m, m);
U = zeros(m, n);
pivcol = [];
tol = sqrt(eps);
sign = 1;

p = 1;
for k = 1:min(m, n)
   [r, p] = findpiv(A(k:m, p:n), k, p, tol);
   if r ~= k
      A([r k], 1:n) = A([k r], 1:n);
      if k > 1, L([r k], 1:k-1) = L([k r], 1:k-1); end
      P([r k], 1:m) = P([k r], 1:m);
      sign = -sign;
   end
   if abs(A(k, p)) >= tol
      pivcol = [pivcol p];
      for i = k+1:m
         L(i, k) = A(i, p) / A(k, p);
         for j = k+1:n
            A(i,j) = A(i, j) - L(i, k)*A(k, j);
         end
      end
   end
   for j = k:n
      U(k, j) = A(k, j) * (abs(A(k, j)) >= tol);
   end
   if p < n, p = p+1; end
end

if nargout < 4
   nopiv = 1:n;
   nopiv(pivcol) = [];
   if ~isempty(pivcol), disp('Pivots in columns:'), disp(pivcol); end
   if ~isempty(nopiv), disp('No pivots in columns:'), disp(nopiv); end
   rank = length(pivcol);
   if rank > 0
      roworder = P*(1:m)';
      disp('Pivots in rows:'), disp(roworder(1:rank)'); end
   end
end
