function [ROW, N, COL, LN] = fourbase(A)

% fourbase  Bases for all 4 fundamental subspaces. 
% 
% [ROW, N, COL, LN] = fourbase(A) returns matrices whose columns
% are bases for the row space, nullspace, column space and
% left nullspace of a matrix A.
%
% The bases for all 4 subspaces come from E*A = R.
% The matrix ROW comes from the r pivot rows of R.
% The matrix N comes from the n-r special solutions. 
% The matrix COL contains the r pivot columns of A.
% The matrix LN contains the last m-r rows of E.
% Those m-r rows multiply A to give the m-r zero rows in R.
% Notice that A = COL * ROW' .
%
% See also rowbasis, nulbasis, colbasis, leftnull.

[m, n] = size(A);
r = rank(A);
E = elim(A);
[R, pivcol] = rref(A);
ROW = R(1:r, :)'; 
N = nulbasis(R);
COL = A(:, pivcol);
LN = E((r+1):m, :)';


