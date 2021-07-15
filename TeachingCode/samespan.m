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
