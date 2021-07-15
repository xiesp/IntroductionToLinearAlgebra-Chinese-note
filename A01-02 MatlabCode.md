# 说明

目前全部来自MIT LA课的Teaching code，后面补充,这里作为注释和学习用，完整代码在 `TachingCode` 目录当中，可用matlab进入目录运行





# 代码



**cab.m**



```matlab
function [c, a, b] = cab(A)

% cab  A = c a b echelon factorization.
%
% [c, a, b] = cab(A) gives echelon bases for the column space in c
% and the row space in b
% b contains the nonzero rows of the echelon form rref(A)
% c contains the nonzero columns of the echelon form rref(A')'
% All extra nonzeros are below I in c and to the right of I in b
% a is the nonsingular submatrix formed by the pivot columns and 
% pivot rows of A.  Those columns of b and rows of c contain I.
%
% See also elim, rref.

[R, pivcol] = rref(A);
[S, pivrow] = rref(A');
b = R(1:rank(A), : );
c = S(1:rank(A), : )';
a = A(pivrow, pivcol);
```





**cafactor.m**

```matlab
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

```



**colbasis.m**

```matlab
function C = colbasis(A)

% colbasis  Basis for the column space. 
%
% C = colbasis(A) returns the r pivot columns of A
% as a basis for the column space of A.
%
% See also fourbase.

[R, pivcol] = rref(A);
C = A(:, pivcol);

```

**cramer.m**

```matlab
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

```

**determ.m**

```matlab
function det = determ(A)

% determ  Matrix determinant from plu.
%
% det = determ(A) computes the determinant of the square matrix A
% as the sign of the permutation times the product of pivots.

[P, L, U, sign] = splu(A);
det = sign * prod(diag(U));

```



**eigen2.m**

```matlab
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

```





**eigshow.m**

```matlab
function eigshow(arg)
%EIGSHOW Graphical demonstration of eigenvalues and singular values.
%
%   EIGSHOW presents a graphical experiment showing the effect on the
%   the unit circle of the mapping induced by various 2-by-2 matrices.
%   A pushbutton allows the choice of "eig" mode or "svd" mode.
%
%   In eig mode, the mouse can be used to move the vector x around the
%   unit circle.  The resulting trajectory of A*x is plotted.  The object
%   is to find vectors x so that A*x is parallel to x.  Each such x is an
%   eigenvector of A.  The length of A*x is the corresponding eigenvalue.
%
%   In svd mode, the mouse moves two perpendicular unit vectors, x and y.
%   The resulting A*x and A*y are plotted.  When A*x is perpendicular to
%   A*y, then x and y are right singular vectors, A*x and A*y are
%   multiples of left singular vectors, and the lengths of A*x and A*y
%   are the corresponding singular values.
%
%   The figure title is a menu of selected matrices, including some
%   with fewer than two real eigenvectors.  EIGSHOW(A) inserts A,
%   which must be 2-by-2, in the menu.
%
%   Here are some questions to consider:
%      Which matrices are singular?
%      Which matrices have complex eigenvalues?
%      Which matrices have double eigenvalues?
%      Which matrices have eigenvalues equal to singular values?
%      Which matrices have nondiagonal Jordan canonical forms?

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1997/11/21 23:25:37 $

if nargin == 0;
initialize
elseif arg == 0
action
elseif arg < 0
setmode(arg)
else
initialize(arg);
end

%------------------

function initialize(arg)

if nargin == 0
arg = 6;
end

if isequal(get(gcf,'tag'),'eigshow');
h = get(gcf,'userdata');
mats = h.mats;
else
set(gcf,'numbertitle','off','menubar','none')
h.svd = 0;
mats = {
'[5/4 0; 0 3/4]'
'[5/4 0; 0 -3/4]'
'[1 0; 0 1]'
'[0 1; 1 0]'
'[0 1; -1 0]'
'[1 3; 4 2]/4'
'[1 3; 2 4]/4'
'[3 1; 4 2]/4'
'[3 1; -2 4]/4'
'[2 4; 2 4]/4'
'[2 4; -1 -2]/4'
'[6 4; -1 2]/4'
'randn(2,2)'};
end

if all(size(arg)==1)
if (arg < length(mats))
mindex = arg;
A = eval(mats{mindex});
else
A = randn(2,2);
S = ['[' sprintf('%4.2f %4.2f; %4.2f %4.2f',A) ']'];
mindex = length(mats);
mats = [mats(1:mindex-1); {S}; mats(mindex)];
end
else
A = arg;
if isstr(A)
S = A;
A = eval(A);
else
S = ['[' sprintf('%4.2f %4.2f; %4.2f %4.2f',A) ']'];
end
if any(size(A) ~= 2)
error('Matrix must be 2-by-2')
end
mats = [{S};  mats];
mindex = 1;
end

clf
if h.svd, t = 'svd / (eig)';
else, t = 'eig / (svd)';
end
uicontrol(...
'style','pushbutton', ...
'units','normalized', ...
'position',[.86 .60 .12 .06], ...
'string',t, ...
'value',h.svd, ...
'callback','eigshow(-1)');
uicontrol(...
'style','pushbutton', ...
'units','normalized', ...
'position',[.86 .50 .12 .06], ...
'string','help', ...
'callback','helpwin eigshow')
uicontrol(...
'style','pushbutton', ...
'units','normalized', ...
'position',[.86 .40 .12 .06], ...
'string','close', ...
'callback','close(gcf)')
uicontrol(...
'style','popup', ...
'units','normalized', ...
'position',[.28 .92 .48 .08], ...
'string',mats, ...
'tag','mats', ...
'fontname','courier', ...
'fontweight','bold', ...
'fontsize',14, ...
'value',mindex, ...
'callback','eigshow(get(gco,''value''))');

s = 1.1*max(1,norm(A));
axis([-s s -s s])
axis square
xcolor = [0 .6 0];
Axcolor = [0 0 .8];
h.A = A;
h.mats = mats;
h.x = initv([1 0]','x',xcolor);
h.Ax = initv(A(:,1),'Ax',Axcolor);
if h.svd
h.y = initv([0 1]','y',xcolor);
h.Ay = initv(A(:,2),'Ay',Axcolor);
xlabel('Make A*x perpendicular to A*y','fontweight','bold')
set(gcf,'name','svdshow')
else
xlabel('Make A*x parallel to x','fontweight','bold')
set(gcf,'name','eigshow')
end
set(gcf,'tag','eigshow', ...
'userdata',h, ...
'windowbuttondownfcn', ...
'eigshow(0); set(gcf,''windowbuttonmotionfcn'',''eigshow(0)'')', ...
'windowbuttonupfcn', ...
'set(gcf,''windowbuttonmotionfcn'','''')')

%------------------

function h = initv(v,t,color)
h.mark = line(v(1),v(2),'marker','.','erase','none','color',color);
h.line = line([0 v(1)],[0 v(2)],'erase','xor','color',color);
h.text = text(v(1)/2,v(2)/2,t,'fontsize',12,'erase','xor','color',color);

%------------------

function action
h = get(gcf,'userdata');
pt = get(gca,'currentpoint');
x = pt(1,1:2)';
x = x/norm(x);
movev(h.x,x);
A = h.A;
movev(h.Ax,A*x);
if h.svd
y = [-x(2); x(1)];
movev(h.y,y);
movev(h.Ay,A*y);
end

%------------------

function movev(h,v)
set(h.mark,'xdata',v(1),'ydata',v(2));
set(h.line,'xdata',[0 v(1)],'ydata',[0 v(2)]);
set(h.text,'pos',v/2);

%------------------

function setmode(arg)
h = get(gcf,'userdata');
h.svd = ~h.svd;
set(gcf,'userdata',h)
initialize(get(findobj(gcf,'tag','mats'),'value'))



```









## 01-02

### findpiv.m

```matlab
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

```



### **elim.m**

```matlab
function [E, R] = elim(A)

% elim  E*A = R factorization.
%
% E = elim(A) returns the elimination matrix E
% that gives the reduced row echelon form E*A = R.
% If A is square and invertible, then E = inv(A).
%
% [E, R] = elim(A) returns the elimination matrix E 
% and the reduced row echelon form R.
%
% See also lu, slu, splu, plu.

[m, n] = size(A);
I = eye(m);
%
% Elimination on the augmented matrix [A I] yields [R E].
%
RE = rref([A I]);
R = RE(:, 1:n);
E = RE(:, (n+1):(m+n));

```



### slu.m

```matlab
function [L, U] = slu(A)

% slu  LU factorization of a square matrix using *no row exchanges*.
%
% [L, U] = slu(A) uses Gaussian elimination to compute a unit 
% lower triangular L and an upper triangular U so that L*U = A.
% The algorithm will stop if a pivot entry is very small.
%
% See also slv, plu, lu.

[n, n] = size(A);

for k = 1:n
if abs(A(k, k)) < sqrt(eps)
disp(['Small pivot encountered in column ' int2str(k) '.'])
end
L(k, k) = 1;
for i = k+1:n
L(i,k) = A(i, k) / A(k, k);
for j = k+1:n
A(i, j) = A(i, j) - L(i, k)*A(k, j);
end
end
for j = k:n
U(k, j) = A(k, j);
end
end

```



### slv.m

```matlab
function x = slv(A, b)
%
% x = slv(A, b)
%
% x = slv(A, b) tries to use the LU factorization computed 
% by slu(A) to solve the linear equation A*x = b.
% Since slu does *no row exchanges*, slv may fail if a 
% small pivot is encountered.
%
% See also slu, solve, slash, partic.

[L, U] = slu(A);

% Forward elimination to solve L*c = b.
% L is lower triangular with 1's on the diagonal.

[n, n] = size(A);
c = zeros(n, 1);
for k = 1:n
s = 0;
for j = 1:k-1
s = s + L(k, j)*c(j);
end
c(k) = b(k) - s;
end

% Back substitution to solve U*x = c.
% U is upper triangular with nonzeros on the diagonal.

x = zeros(n, 1);
for k = n:-1:1
t = 0;
for j = k+1:n
t = t + U(k, j)*x(j);
end
x(k) = (c(k) - t) / U(k, k);
end
x = x';

```





### splu.m

```matlab
function [P, L, U, sign] = splu(A)

% splu  Square PA=LU factorization *with row exchanges*.
%
% [P, L, U] = splu(A), for a square, invertible matrix A,
% uses Gaussian elimination to compute a permutation
% matrix P, a lower triangular matrix L and 
% an upper triangular matrix U so that P*A = L*U.
% P, L and U are the same size as A.
% sign = det(P); it is 1 or -1.
%
% See also slu, lu, rref, partic, nulbasis, colbasis.

[m, n] = size(A);
if m ~= n
error('Matrix must be square.')
end
P = eye(n, n);
L = eye(n, n);
U = zeros(n, n);
tol = sqrt(eps);
sign = 1;

for k = 1:n
if abs(A(k, k)) < tol
for r = k:n
if abs(A(r, k)) >= tol
break
end
if r == n
if nargout == 4
sign = 0;
return
else
disp('A is singular within tolerance.')
error(['No pivot in column ' int2str(k) '.'])
end
end
end
A([r k], 1:n) = A([k r], 1:n);
if k > 1, L([r k], 1:k-1) = L([k r], 1:k-1); end
P([r k], 1:n) = P([k r], 1:n);
sign = -sign;
end
for i = k+1:n
L(i, k) = A(i, k) / A(k, k);
for j = k+1:n
A(i, j) = A(i, j) - L(i, k)*A(k, j);
end
end
for j = k:n
U(k, j) = A(k, j) * (abs(A(k, j)) >= tol);
end
end

if nargout < 4
roworder = P*(1:n)';
disp('Pivots in rows:'), disp(roworder'); end
end

```





## 01-03


### samespan.m

```matlab
function samespan(A1, A2)

% samespan  Test if two matrices have the same column space.
%
% samespan(A1, A2) 
% If the column spaces of A1 and A2 are the same,
% the function returns the dimension of this subspace.
% If the subspaces are different, the function returns 0.
%
% See also rank.

rankA1 = rank(A1)
rankA2 = rank(A2)
rankboth = rank([A1 A2])

if (rankA1 == rankA2) & (rankA1 == rankboth) 
disp('A1 and A2 have the same column space.');
else
disp('A1 and A2 have different column spaces.');
end

```





**fourbas.m**

```matlab
function [ROW, N, COL, LN] = fourbase(A)

% fourbase  Bases for all 4 fundamental subspaces. 
% 
% [ROW, N, COL, LN] = fourbase(A) returns matrices whose columns
% are bases for the row space, nullspace, column space and
% left nullspace of a matrix A.
%
% The bases for all 4 subspaces come from E*A = R.
% The matrix ROW comes from the r pivot rows of R.
% The matrix N comes from the n-r special solutions. 
% The matrix COL contains the r pivot columns of A.
% The matrix LN contains the last m-r rows of E.
% Those m-r rows multiply A to give the m-r zero rows in R.
% Notice that A = COL * ROW' .
%
% See also rowbasis, nulbasis, colbasis, leftnull.

[m, n] = size(A);
r = rank(A);
E = elim(A);
[R, pivcol] = rref(A);
ROW = R(1:r, :)'; 
N = nulbasis(R);
COL = A(:, pivcol);
LN = E((r+1):m, :)';



```





## 01-04

### grams.m

```matlab
function [Q, R] = grams(A)

% grams  Gram-Schmidt orthogonalization of the columns of A.
% The columns of A are assumed to be linearly independent.
%
% Q = grams(A) returns an m by n matrix Q whose columns are 
% an orthonormal basis for the column space of A.
%
% [Q, R] = grams(A) returns a matrix Q with orthonormal columns
% and an invertible upper triangular matrix R so that A = Q*R.
%
% Warning: For a more stable algorithm, use [Q, R] = qr(A, 0) .

[m, n] = size(A);
Asave = A;
for j = 1:n
for k = 1:j-1
mult = (A(:, j)'*A(:, k)) / (A(:, k)'*A(:, k));
A(:, j) = A(:, j) - mult*A(:, k);
end
end
for j = 1:n
if norm(A(:, j)) < sqrt(eps)
error('Columns of A are linearly dependent.')
end
Q(:, j) = A(:, j) / norm(A(:, j));
end
R = Q'*Asave;

```



**house.m**

```
function X = house

% X = house stores the "house" data set in X.
%
% The "house" data set is for use with plot2d.
% Try plot2d(A*X) for various 2 by 2 matrices A.
%
% See also hand.

X = [ -6  -6  -7   0   7   6   6  -3  -3   0   0  -6
-7   2   1   8   1   2  -7  -7  -2  -2  -7  -7 ];

```



**inverse.m**

```matlab
function Ainv = inverse(A)

% inverse  Matrix inverse by Gauss-Jordan elimination.
%
% Ainv = inverse(A) computes the inverse of A, if it exists.
%
% Row reduction applied to [A I] using elim produces [I Ainv].
%
% See also inv, elim. 

[m, n] = size(A);
r = rank(A);
if (r == m) & (r == n) 
[Ainv, R] = elim(A);
else
Ainv = [];
disp('Warning: A is not a square, invertible matrix.');
end;

```



**leftnull.m**

```matlab
function LN = leftnull(A)

% leftnull  Basis for the left nullspace.
%
% LN = leftnull(A) returns a basis for the 
% left nullspace in the *columns* of LN.
%
% The left nullspace of A is the nullspace of A'.
% The command fourbase(A) finds a different basis
% for the left nullspace of A. 
%
% See also fourbase.

LN = nulbasis(A');

```

**linefit.m**

```matlab
function linefit(t, b)

% linefit  Plot the least squares fit by a line.
%
% linefit(t, b), where t and b are vectors of the same length,
% displays the best line fit to the data points (t(i), b(i)).

% We now insure that t and b are column vectors.
t = t(:); b = b(:);

% Form the matrix whose first column is all ones and
% whose second column is the vector t.
n = length(t);
e = ones(n, 1);
A = [e t];

% Solve the least squares problem, A*x ~= b.
xhat = lsq(A, b);
c = xhat(1);
d = xhat(2);

% Plot the results.
tline = [1.1*min(t)-0.1*max(t), 1.1*max(t)-0.1*min(t)];
yline = c + d*tline;
plot(t,b,'ro',t,c+d*t,'k*',tline,yline,'k-')
if d >= 0 , sign = ' + '; else, sign = ' - '; end
title(['Best line is ' num2str(c) sign num2str(abs(d)) '*t.'])
xlabel('t')

```

**lsq.m**

```matlab
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

```



**normal.m**

```matlab
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

```

**nullbasis.m**

```matlab
function N = nulbasis(A)

% nulbasis  Basis for nullspace.
%
% N = nulbasis(A) returns a basis for the nullspace of A
% in the columns of N. The basis contains the n-r special 
% solutions to Ax=0.  freecol is the list of free columns.
%
% Example:
%
% >> A = [1 2 0 3;
%        [0 0 1 4];
%
% >> N = nulbasis(A)
%
%    N = [-2  -3]   
%        [ 1   0]
%        [ 0  -4]
%        [ 0   1]
%
% See also fourbase.

[R, pivcol] = rref(A, sqrt(eps));
[m, n] = size(A);
r = length(pivcol);
freecol = 1:n;
freecol(pivcol) = [];
N = zeros(n, n-r);
N(freecol, : ) = eye(n-r);
N(pivcol,  : ) = -R(1:r, freecol);

```







**orthcomp.m**

```matlab
function BCOMP = orthcomp(B)

% orthcomp  Orthogonal complement of a subspace.
%
% BCOMP = orthcomp(B) returns a basis for the 
% orthogonal complement of the column space of B.
% This subspace contains all vectors orthogonal
% to the column space of B.
% It is the left nullspace of B.
%
% See also leftnull, nulbasis.

BCOMP = leftnull(B);

```



**partic.m**

```matlab
function x = partic(A, b)

% partic  Particular solution of Ax=b.
%
% x = partic(A, b) returns a particular solution to Ax=b.
% This particular solution has all free variables set to zero.
% An empty vector is returned if Ax=b is not solvable.
%
% See also slash as in A\b .

[m, n] = size(A);
[Rd, pivcol] = rref([A b]);
r = length(pivcol);
%
% If the last column of the augmented matrix [A b] 
% is a pivot column, then Ax=b has no solution.
%
if max(pivcol) == n+1
x = [];
else
%
% The values of the pivot variables are in the
% last column (which is called d) of Rd.
% The free variables are zero in this particular solution.
%
x = zeros(n, 1);
d = Rd(:, n+1);
x(pivcol) = d(1:r);
end

```





**plot2d.m**

```matlab
function plot2d(X)

% plot2d  Two dimensional plot.

% X is a matrix with 2 rows and any number of columns.
% plot2d(X) plots these columns as points in the plane
% and connects them, in order, with lines.
% The scale is set to [-10, 10] in both directions.
%
% For example, the statement X=house creates 
% a sample matrix X representing common figures.
% Then, for various 2 by 2 matrices A,
%    plot2d(A*X)
% demonstrates the effect of multiplication by A.

x = X(1,:)';
y = X(2,:)';
plot(x, y, 'ro', x, y, 'g-');
axis([-10 10 -10 10])
axis('square')

```







**plu.m**

```matlab
function [P, L, U, pivcol, sign] = plu(A)

% plu  Rectangular PA=LU factorization *with row exchanges*.
%
% [P, L, U] = plu(A), for a rectangular matrix A, uses Gaussian elimination
% to compute a permutation matrix P, a lower triangular matrix L and 
% an upper trapezoidal matrix U so that PA = LU.
% U is the same size as A.  
% P and L are square, with as many rows as A.
% sign = det(P); it is 1 or -1.
%
% See also elim, slu, lu, rref, partic, nulbasis, colbasis.

[m, n] = size(A);
P = eye(m, m);
L = eye(m, m);
U = zeros(m, n);
pivcol = [];
tol = sqrt(eps);
sign = 1;

p = 1;
for k = 1:min(m, n)
[r, p] = findpiv(A(k:m, p:n), k, p, tol);
if r ~= k
A([r k], 1:n) = A([k r], 1:n);
if k > 1, L([r k], 1:k-1) = L([k r], 1:k-1); end
P([r k], 1:m) = P([k r], 1:m);
sign = -sign;
end
if abs(A(k, p)) >= tol
pivcol = [pivcol p];
for i = k+1:m
L(i, k) = A(i, p) / A(k, p);
for j = k+1:n
A(i,j) = A(i, j) - L(i, k)*A(k, j);
end
end
end
for j = k:n
U(k, j) = A(k, j) * (abs(A(k, j)) >= tol);
end
if p < n, p = p+1; end
end

if nargout < 4
nopiv = 1:n;
nopiv(pivcol) = [];
if ~isempty(pivcol), disp('Pivots in columns:'), disp(pivcol); end
if ~isempty(nopiv), disp('No pivots in columns:'), disp(nopiv); end
rank = length(pivcol);
if rank > 0
roworder = P*(1:m)';
disp('Pivots in rows:'), disp(roworder(1:rank)'); end
end
end

```



**poly2str.m**

```matlab
function p = poly2str(c, x)

% poly2str  Convert a polynomial coefficient vector to a string.
%
% p = poly2str(c) generates a string representation of the polynomial
% whose coefficents are in the vector c.  
% The default variable is 'x', unless otherwise specified by 
% p = poly2str(c, 's').
% The coefficients are approximated, if necessary, by the rational
% values obtained from rat.
%	
% If x has a numeric value and the elements of c are reproduced
% exactly by rat, then eval(poly2str(c)) will return the same value 
% as polyval(c, x).
%
% See also polyval, rat.

if nargin < 2, x = 'x'; end
if all(c == 0), p = '0'; return, end

p = [];
n = length(c);
for d = 0: n-1
if d > 0
if c(n-d+1) > 0
p = [' + ' p];
elseif c(n-d+1) < 0
p = [' - ' p];
end
end
if c(n-d) ~= 0
if d == 1
p = [x p];
elseif d > 1
p = [x '^' int2str(d) p];
end
if (abs(c(n-d)) ~= 1) | (d==0)
if d > 0,
p = ['*' p];
end
[sn, sd] = rat(abs(c(n-d)));
s = num2str(sn);
if sd ~= 1, s = [s '/' num2str(sd)]; end
p = [s p];
end
end
end
if n > 0
if c(1) < 0
p = ['-' p];
end
end

```



**project.m**

```matlab
function [p, e] = project(A, b)

% project  Project a vector b onto the column space of A.
%
% p = project(A, b) returns the orthogonal projection of a 
% vector b onto the column space of A.
%
% [p, e] = project(A, b) also returns the vectors e = b - p.
% p is the projection of b onto the column space of A.
% e is the projection of b onto the left nullspace of A.
% Notice that b = p + e and p' * e = 0. 
%
% See also projmat.

P = projmat(A);
p = P * b;
e = b - p;

```



**projmat.m**

```matlab
function P = projmat(A)

% projmat  Projection matrix for the column space of A.
%
% P = projmat(A) returns the projection matrix for
% the column space of A.

A = colbasis(A);
P = A*inv(A'*A)*A';

```





**randperm.m**

```matlab
function p = randperm(n)

% randperm  Random permutation.
%
% p = randperm(n) returns a random permutation of 1:n.

[ignore, p] = sort(rand(1, n));

```





#### rowbasis.m

```matlab
function B = rowbasis(A)

% rowbasis  Basis for the row space. 
%
% B = rowbasis(A) returns a basis for the row space of A
% in the *columns* of B.
% The row space of A is the column space of A'.
% rowbasis finds the first r linearly independent 
% columns of A'.
%
% The command fourbase(A) uses rref(A) to find a 
% different basis for the row space of A.
%
% See also fourbase.

B = colbasis(A');

```

## 01-06

### eigval.m

```matlab
function [evalues, repeats] = eigval(A)

% eigval  Eigenvalues and their algebraic multiplicity.
%
% evalues = eigval(A) returns the distinct eigenvalues of A,
% with duplicates removed and sorted in decreasing order.
%
% [evalues, repeats] = eigval(A) also returns the row vector
% repeats that gives the multiplicity of each eigenvalue.
% The sum of the multiplicities is n.
%
% Examples: Let A = eye(n) and B = diag([3 4]).
% For A, evalues is 1 and repeats is n.
% For B, evalues is [4; 3]  and repeats is [1 1].

tol = sqrt(eps);
lambda = sort(eig(A));
lambda = round(lambda/tol) * tol;
%
% lambda gives all n eigenvalues (repetitions included).
%
evalues = unique(lambda);
evalues = flipud(evalues);
n = length(lambda);
d = length(evalues);
A = ones(n, 1) * evalues';
B = lambda * ones(1, d);
MATCH = abs(A-B) <= tol;
%
% MATCH is an n by d zero matrix except
% MATCH(i,j) = 1 when lambda(i) = evalues(j).
% Summing the columns gives the row vector repeats.
%
repeats = sum(MATCH);

```



### eigvec.m

```matlab
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

```







**signperm.m**

```matlab
function [sign, PERM] = signperm(p)

% signperm  Determinant of the permutation matrix with rows ordered by p.
%
% sign = signperm(p) returns the sign of the 
% permutation associated with the vector p.
%
% [sign, PERM] also returns the permutation matrix PERM.
%
% Example: Let p = [2 3 1].
% Then sign = 1 and PERM = [0 1 0; 0 0 1; 1 0 0] .
% 

n = length(p);
I = eye(n);
PERM = I(p, :);
sign = det(PERM);

```







#### splv.m

```matlab
function x = splv(A, b)

% splv  The solution to a square, invertible system.
% x = splv(A, b) uses the PA = LU factorization
% computed by splu to solve Ax = b.
%
% See also slv, splu, slash.

[P, L, U] = splu(A);
[n, n] = size(A);

% Permute the right hand side.
b = P*b;

% Forward elimination to solve L*c = b.
c = zeros(n, 1);
for k = 1:n
s = 0;
for j = 1:k-1
s = s + L(k, j)*c(j);
end
c(k) = b(k) - s;
end

% Back substitution to solve U*x = c.
x = zeros(n, 1);
for k = n:-1:1
t = 0;
for j = k+1:n
t = t + U(k, j)*x(j);
end
x(k) = (c(k) - t) / U(k, k);
end

```



#### symmeig.m

```matlab
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

```



#### tridiag.m

```matlab
function T = tridiag(a, b, c, n)

% tridiag  Tridiagonal matrix.
% T = tridiag(a, b, c, n) returns an n by n matrix that has 
% a, b, c as the subdiagonal, main diagonal, and superdiagonal 
% entries in the matrix.

T = b*diag(ones(n,1)) + c*diag(ones(n-1,1),1) + a*diag(ones(n-1,1),-1);

```
