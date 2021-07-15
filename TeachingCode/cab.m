function [c, a, b] = cab(A)

% cab  A = c a b echelon factorization.
%
% [c, a, b] = cab(A) gives echelon bases for the column space in c
% and the row space in b
% b contains the nonzero rows of the echelon form rref(A)
% c contains the nonzero columns of the echelon form rref(A')'
% All extra nonzeros are below I in c and to the right of I in b
% a is the nonsingular submatrix formed by the pivot columns and 
% pivot rows of A.  Those columns of b and rows of c contain I.
%
% See also elim, rref.

[R, pivcol] = rref(A);
[S, pivrow] = rref(A');
b = R(1:rank(A), : );
c = S(1:rank(A), : )';
a = A(pivrow, pivcol);
