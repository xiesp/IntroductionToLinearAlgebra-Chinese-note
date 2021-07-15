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
