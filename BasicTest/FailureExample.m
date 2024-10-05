% This script vaildates the failure transition matrix of poisson failure model

close all
clear all
clc
% rng('default')

% Test Param
dt_sim = 1; % day
tf = 100; % day
time_sim = 1:dt_sim:tf; % 365 day
lambda = 15; % # failure per year
N_sat = 40; % Initial # of sat.
iter_max = 10000; % # of iteration

% Run simulation
lambda_sim  = lambda/365; % failure per day
N_fail = CustomPoisRnd(lambda_sim, length(time_sim), iter_max); % sample # of failure
N_on = N_sat - cumsum(N_fail,1); % # of sat. at every day 
N_on(N_on < 0) = 0;

% Analysis
nx = N_sat + 1; % Dim. of state
Pf = zeros(nx); % Failure transition matrix
for i = 1:nx
    Pf(i:end,i) = poisspdf(0:nx-i, lambda_sim)';
    Pf(end,i) = 1 - sum(Pf(1:end-1, i));
end
pi_0 = zeros(nx,1);
pi_0(1) = 1;
pi_f = Pf^tf*pi_0;

% Result Plot
figure(1); hold on % Time series graph
plot(time_sim, N_on)
plot(time_sim, mean(N_on,2),'k-','linewidth', 2)
xlabel('Day'); ylabel('# of sat.')

figure(2); hold on % Distribution at the last day
histogram(N_on(end,:),'Normalization','probability')
plot((N_sat:-1:0), pi_f', 'r*')


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