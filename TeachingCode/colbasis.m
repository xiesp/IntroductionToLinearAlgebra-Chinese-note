function C = colbasis(A)

% colbasis  Basis for the column space. 
%
% C = colbasis(A) returns the r pivot columns of A
% as a basis for the column space of A.
%
% See also fourbase.

[R, pivcol] = rref(A);
C = A(:, pivcol);
