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
