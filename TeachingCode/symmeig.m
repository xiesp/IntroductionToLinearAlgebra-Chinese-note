function [Q, LAMBDA] = symmeig(A)

% symmeig  Eigenvalues and eigenvectors of a symmetric matrix.
% The matrix A is assumed to be symmetric.
%
% Q = symmeig(A) returns a set of orthonormal eigenvectors for A.
%
% [Q, LAMBDA] = symmeig(A) also returns the corresponding eigenvalues 
% on the diagonal of LAMBDA. The eigenvalues in LAMBDA are in 
% decreasing order.
%
% See also eigval, eigvec, eig.

[m, n] = size(A);
if norm(A'-A) <= sqrt(eps)
%
% The eigenvectors in S are linearly independent but not orthogonal.
% Eigenvectors for different eigenvalues *are* orthogonal.
% Gram-Schmidt (qr(S)) gives orthonormal eigenvectors in Q.
%
  [S, LAMBDA] = eigvec(A);
  [Q, R] = qr(S);
else
  Q = []; LAMBDA = [];
  error('The matrix is not symmetric.');
end
