function [xhat, p, e] = lsq(A, b)

% lsq  Least squares solution to Ax=b.
% The columns of A are linearly independent.
%
% [xhat, p, e] = lsq(A, b) finds a least squares
% solution to the overdetermined system Ax ~= b.
% xhat solves the normal equations A'*A*xhat = A'*b.
% p is the projection of b onto the column space.
% e = b - p.

if det(A'*A) == 0
   error('Columns of A are linearly dependent.')
end
xhat = partic(A'*A, A'*b);
p = A * xhat;
e = b - p;
