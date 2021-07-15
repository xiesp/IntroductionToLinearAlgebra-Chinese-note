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
