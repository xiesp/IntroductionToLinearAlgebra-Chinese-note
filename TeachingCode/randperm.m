function p = randperm(n)

% randperm  Random permutation.
%
% p = randperm(n) returns a random permutation of 1:n.

[ignore, p] = sort(rand(1, n));
