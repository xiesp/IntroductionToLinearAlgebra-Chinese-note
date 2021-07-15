function T = tridiag(a, b, c, n)

% tridiag  Tridiagonal matrix.
% T = tridiag(a, b, c, n) returns an n by n matrix that has 
% a, b, c as the subdiagonal, main diagonal, and superdiagonal 
% entries in the matrix.

T = b*diag(ones(n,1)) + c*diag(ones(n-1,1),1) + a*diag(ones(n-1,1),-1);
