% This script tests the custom-made poisson random number generator
% To validate the custom-fucntion, you will need the toolbox for poissrnd function

clear all
close all
clc

lamda = 2; % Dist. Parameter
n = 1000; % Number of function call

tic
for i = 1:n
    Zmy = CustomPoisRnd(lamda,10000,1);
end
Tmy = toc;

tic
for i = 1:n
    Zmat = poissrnd(lamda,10000,1);
end
Tmat = toc;

disp([Tmy, Tmat]/n)

edges = 0:1:11;
h1 = histcounts(Zmy,edges);
h2 = histcounts(Zmat,edges);

figure(1)
bar(edges(1:end-1),[h1; h2])


function [Y] = CustomPoisRnd(p, m, n)
    % Input Handling
    if nargin == 1
        m = 1;
        n = 1;
    elseif nargin == 2
        n = 1;
    end
    
    % Compute Poisson CDF (consider upto mean + 10-sigma)
    Zset = 0:1:ceil(11*p);
    PDF = p.^Zset*exp(-p)./factorial(Zset); 
    CDF = cumsum(PDF);
    CDF(end) = 1;
    
    % Generate uniform distribution
    Z = rand(m,n);
    Y = zeros(m,n);
    for i = 1:m
        for j = 1:n
            tmp = CDF - Z(i,j);
            idx = find(tmp>=0,1);
            Y(i,j) = Zset(idx);
        end
    end
end