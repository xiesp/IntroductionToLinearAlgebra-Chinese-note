function Ainv = inverse(A)

% inverse  Matrix inverse by Gauss-Jordan elimination.
%
% Ainv = inverse(A) computes the inverse of A, if it exists.
%
% Row reduction applied to [A I] using elim produces [I Ainv].
%
% See also inv, elim. 

[m, n] = size(A);
r = rank(A);
if (r == m) & (r == n) 
  [Ainv, R] = elim(A);
else
  Ainv = [];
  disp('Warning: A is not a square, invertible matrix.');
end;
