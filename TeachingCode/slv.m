function x = slv(A, b)
%
% x = slv(A, b)
%
% x = slv(A, b) tries to use the LU factorization computed 
% by slu(A) to solve the linear equation A*x = b.
% Since slu does *no row exchanges*, slv may fail if a 
% small pivot is encountered.
%
% See also slu, solve, slash, partic.

[L, U] = slu(A);

% Forward elimination to solve L*c = b.
% L is lower triangular with 1's on the diagonal.

[n, n] = size(A);
c = zeros(n, 1);
for k = 1:n
   s = 0;
   for j = 1:k-1
      s = s + L(k, j)*c(j);
   end
   c(k) = b(k) - s;
end

% Back substitution to solve U*x = c.
% U is upper triangular with nonzeros on the diagonal.

x = zeros(n, 1);
for k = n:-1:1
   t = 0;
   for j = k+1:n
      t = t + U(k, j)*x(j);
   end
   x(k) = (c(k) - t) / U(k, k);
end
x = x';
