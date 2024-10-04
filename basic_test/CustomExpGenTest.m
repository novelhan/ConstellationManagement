% This script tests the custom-made exponential random number generator
% To validate the custom-fucntion, you will need the toolbox for exprnd function

clear all
close all
clc

mu = 2; % Dist. Parameter
n = 1000; % Number of function call

tic
for i = 1:n
    Zmy = CustomExpRnd(mu,10000,1);
end
Tmy = toc;

tic
for i = 1:n
    Zmat = exprnd(mu,10000,1);
end
Tmat = toc;

disp([Tmy, Tmat]/n)

edges = 0:1:15;
h1 = histcounts(Zmy,edges);
h2 = histcounts(Zmat,edges);

figure(1)
bar(edges(1:end-1),[h1; h2])


function [Y] = CustomExpRnd(mu, m, n)
    % Input Handling
    if nargin == 1
        m = 1;
        n = 1;
    elseif nargin == 2
        n = 1;
    end
    
    % Mean = 1/lambda
    p = 1/mu;
    
    % Generate uniform distribution
    Z = rand(m,n);
    
    % Apply Inverse Transform
    Y = -1/p*log(1-Z);
end