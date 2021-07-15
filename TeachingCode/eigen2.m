function eigen2(A)

% eigen2  Characteristic polynomial, eigenvalues, eigenvectors 
% of a 2 by 2 matrix.
%
% eigen2(A) prints the characteristic polynomial det(A-e*I),
% eigenvalues, and eigenvectors of A. 
%
% If A is not diagonalizable, its single eigenvector is 
% printed twice.

d = A(1,1)*A(2,2) - A(1,2)*A(2,1);
t = A(1,1) + A(2,2);
e1 = (t + sqrt(t^2 - 4*d))/2;
e2 = (t - sqrt(t^2 - 4*d))/2;
if A(1,2) ~= 0
   x1 = [A(1,2); e1-A(1,1)];
   x2 = [A(1,2); e2-A(1,1)];
elseif A(2,1) ~= 0
   x1 = [e1-A(2,2); A(2,1)];
   x2 = [e2-A(2,2); A(2,1)];
else
   x1 = [1; 0];
   x2 = [0; 1];
end

disp(' ')
disp('For this matrix, the polynomial whose roots are the eigenvalues is:')
disp(['   e^2 - ' num2str(t) '*e + ' num2str(d) ' = 0'])

disp(' ')
disp('The first eigenvalue and eigenvector are:')
e1
x1

disp(' ')
disp('The second eigenvalue and eigenvector are:')
e2
x2
