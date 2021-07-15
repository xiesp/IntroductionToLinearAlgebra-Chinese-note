function x = partic(A, b)

% partic  Particular solution of Ax=b.
%
% x = partic(A, b) returns a particular solution to Ax=b.
% This particular solution has all free variables set to zero.
% An empty vector is returned if Ax=b is not solvable.
%
% See also slash as in A\b .

[m, n] = size(A);
[Rd, pivcol] = rref([A b]);
r = length(pivcol);
%
% If the last column of the augmented matrix [A b] 
% is a pivot column, then Ax=b has no solution.
%
if max(pivcol) == n+1
  x = [];
else
%
% The values of the pivot variables are in the
% last column (which is called d) of Rd.
% The free variables are zero in this particular solution.
%
  x = zeros(n, 1);
  d = Rd(:, n+1);
  x(pivcol) = d(1:r);
end
