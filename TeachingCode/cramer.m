function x = cramer(A, b)

% cramer  Solve the system Ax=b.
% The matrix A is square and invertible.
%
% x = cramer(A, b) solves the square system Ax = b.

[m, n] = size(A);
if m ~= n
  error('Matrix is not square.') 
end
if det(A) == 0
  error('Matrix is singular.')
end
for j = 1:n
  B = A;
  B(:, j) = b;
  x(j) = det(B) / det(A);
end
x = x';
