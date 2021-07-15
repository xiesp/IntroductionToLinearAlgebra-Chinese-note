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
