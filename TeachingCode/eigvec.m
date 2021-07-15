function [S, D] = eigvec(A)

% eigvec  Eigenvectors and their geometric multiplicity.
%
% S = eigvec(A) returns the largest possible set of linearly
% independent eigenvectors of A. 
%
% [S, D] = eigvec(A) also returns the corresponding eigenvalues
% in the diagonal matrix D.
% Each eigenvalue in D is repeated according to the number of its
% linearly independent eigenvectors. This is its geometric multiplicity.
%
% Always A*S = S*D. If S is square then A is diagonalizable and
% inv(S)*A*S = D = LAMBDA.

[m, n] = size(A);
I = eye(n);
[evalues, repeats] = eigval(A);
S = []; d = [];
for k = 1 : length(evalues);
  s = nulbasis(A - evalues(k)*I);
  [ms, ns] = size(s);
  S = [S s];
  temp = ones(ns, 1) * evalues(k);
  d = [d; temp];
end
D = diag(d);
