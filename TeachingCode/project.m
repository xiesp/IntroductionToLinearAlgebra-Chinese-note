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
