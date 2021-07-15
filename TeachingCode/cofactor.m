function C = cofactor(A, i, j)

% cofactor  Matrix of cofactors.
%
% C = cofactor(A) returns the matrix of cofactors of A.
% If A is invertible, then inv(A) = C' / det(A).
%
% C = cofactor(A, i, j) returns the cofactor of 
% row i, column j of A.

if nargin == 3
% Remove row i and column j to produce the minor.
  M = A;
  M(i,:) = [];
  M(:,j) = [];
  C = (-1)^(i+j) * det(M);
else
  [n,n] = size(A);
  for i = 1:n
    for j = 1:n
    C(i,j) = cofactor(A, i, j);
    end
  end
end
