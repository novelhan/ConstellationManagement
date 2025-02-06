% This function computes the P.D.F of Poisson distribution. It is defined to replace 'poisspdf'
% fucntion which is toolbox dependant.
%
% Input: 
% k: Random variable for # of event (can be a scalar/vector/matrix)
% p: Pois. parameter (# of event per unit length)
%
% Output:
% pdf: P.D.F of Pois. under the inputs
%
% Reference
% https://en.wikipedia.org/wiki/Poisson_distribution

function [pdf] = CustomPoisPdf(k,p)
    pdf = p.^k*exp(-p)./factorial(k);
end

