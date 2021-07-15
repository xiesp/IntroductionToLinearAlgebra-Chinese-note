function [U, LAMBDA] = normal(A)

% normal  Eigenvalues and eigenvectors of a normal matrix A.
%
% U = normal(A) returns a set of orthonormal eigenvectors for A.
%
% [U, LAMBDA] = normal(A) also returns the corresponding eigenvalues 
% on the diagonal of LAMBDA. The eigenvalues on the diagonal of 
% LAMBDA are sorted by magnitude.
%
% Normal matrices (A'*A = A*A') may have complex eigenvalues and 
% eigenvectors. If A itself is complex, A' is its conjugate transpose.
%
% See also eig, eigval, eigvec, symmeig.

[m, n] = size(A);
E = A'*A - A*A';
if norm(E) <= sqrt(eps)
%
% The eigenvectors in S are linearly independent but not orthogonal.
% Eigenvectors for different eigenvalues *are* orthogonal.
% Gram-Schmidt (qr(S)) gives orthonormal eigenvectors in U.
%
  [S, LAMBDA] = eigvec(A);
  [U, R] = qr(S);
else
  U = []; LAMBDA = [];
  error('The matrix is not normal.');
end
