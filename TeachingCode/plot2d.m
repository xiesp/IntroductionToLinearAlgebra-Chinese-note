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
