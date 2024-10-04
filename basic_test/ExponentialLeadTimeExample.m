% This script vaildates the exact probability equation for the shifted exponential distribution

clear all
close all
clc

% rng('default')
% Test Param
dt_sim = 1; % day
mu_LV = 45; % day
b_LV = 30; % day
iter_max = 50000; % # of iteration

% Run simulation
T_LV = b_LV + CustomExpRnd(mu_LV,iter_max,1);
T_LV = ceil(T_LV); % Apply Ceiling function

% Do analysis
day_max = b_LV + mu_LV*5; % Max day for analysis
rho_LV = exp(-( (b_LV:day_max-1) - b_LV)/mu_LV)*(1-exp(-dt_sim/mu_LV));
rho_LV = [zeros(1,b_LV), rho_LV];

% Plot Result
dT_edge = 0.5:1:(day_max+0.5);
figure(1); hold on
histogram(T_LV, dT_edge(1:end-1), 'Normalization','probability')
plot(1:day_max, rho_LV, 'r*')


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