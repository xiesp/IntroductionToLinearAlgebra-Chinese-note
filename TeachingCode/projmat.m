function P = projmat(A)

% projmat  Projection matrix for the column space of A.
%
% P = projmat(A) returns the projection matrix for
% the column space of A.

A = colbasis(A);
P = A*inv(A'*A)*A';
