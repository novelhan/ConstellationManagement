function x = limitdist1(P)
%Obtain the stationary probability distribution
%vector p of an irreducible, recurrent Markov
%chain by state reduction. P is the transition
%probabilities matrix of a discrete-time Markov
%chain or the generator matrix Q.
% https://www.math.wustl.edu/~feres/Math450Lect04.pdf

[n, ~] = size(P);
A = [eye(n)-P; ones(1,n)];
b = [zeros(n,1); 1];
x = A\b;

end