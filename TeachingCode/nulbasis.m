function N = nulbasis(A)

% nulbasis  Basis for nullspace.
%
% N = nulbasis(A) returns a basis for the nullspace of A
% in the columns of N. The basis contains the n-r special 
% solutions to Ax=0.  freecol is the list of free columns.
%
% Example:
%
% >> A = [1 2 0 3;
%        [0 0 1 4];
%
% >> N = nulbasis(A)
%
%    N = [-2  -3]   
%        [ 1   0]
%        [ 0  -4]
%        [ 0   1]
%
% See also fourbase.

[R, pivcol] = rref(A, sqrt(eps));
[m, n] = size(A);
r = length(pivcol);
freecol = 1:n;
freecol(pivcol) = [];
N = zeros(n, n-r);
N(freecol, : ) = eye(n-r);
N(pivcol,  : ) = -R(1:r, freecol);
