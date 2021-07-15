function [k, p] = findpiv(A, k, p, tol)

% findpiv  Used by plu to find a pivot for Gaussian elimination.
%
% [r, p] = findpiv(A(k:m, p:n), k, p, tol) finds the first element in
% the specified submatrix which is larger than tol in absolute value.
% It returns indices r and p so that A(r, p) is the pivot.

[m,  n] = size(A);
r = find(abs(A(:)) > tol);
if isempty(r), return, end
%
r = r(1);
j = fix((r-1)/m)+1;
p = p+j-1;
k = k+r-(j-1)*m-1;
