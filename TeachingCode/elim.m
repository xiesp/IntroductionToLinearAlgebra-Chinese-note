function [E, R] = elim(A)

% elim  E*A = R factorization.
%
% E = elim(A) returns the elimination matrix E
% that gives the reduced row echelon form E*A = R.
% If A is square and invertible, then E = inv(A).
%
% [E, R] = elim(A) returns the elimination matrix E 
% and the reduced row echelon form R.
%
% See also lu, slu, splu, plu.

[m, n] = size(A);
I = eye(m);
%
% Elimination on the augmented matrix [A I] yields [R E].
%
RE = rref([A I]);
R = RE(:, 1:n);
E = RE(:, (n+1):(m+n));
