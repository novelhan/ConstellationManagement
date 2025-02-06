% This script tests effect of different step-size for finite state Markov system
close all
clear all
clc

%.. PATH
mfilepath = pwd;
idcs = strfind(mfilepath,'\');
libdir = mfilepath(1:idcs(end));
addpath([libdir, 'CommonSource'])

% rng('default')

% Test Param
dt_sim = 1; % day
tf = 30; % day
time_sim = 0:dt_sim:tf; % 365 day
p_day = 15/365; % # failure per day
N0 = 40; % Initial # of sat.
iter_max = 1000; % # of iteration

% Run simulation
p_sim  = p_day*dt_sim; % failure per day

%% Run Each Simulation
% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(N0+1, iter_max);

for itr = 1:iter_max
    % The number of satellite at current time step
    Non_k = N0;       
    
    % Apply Policy
    for k = 1:length(time_sim)
        % Save the state
        Non(k,itr) = Non_k;
        
        % State-dependant failure
        N_fail = CustomPoisRnd(Non_k*p_sim, 1);
        
        % Update the state
        Non_k = max([Non_k - N_fail, 0]);
    end
end

%% Check with various dt_mc
dt_mc = 1; % It seems that if dt_mc is not too different from dt_sim then it gives similar results
p_mc = p_sim*dt_mc;
t_mc = 0:dt_mc:tf;
cf = length(t_mc);
nx = N0+1;
I = eye(nx);
Pf = zeros(nx);
for i = 1:nx
    Pf(i:end,i) = CustomPoisPdf(0:nx-i, (nx-i)*p_mc)';  
    Pf(end,i) = 1 - sum(Pf(1:end-1, i));
end
pi_0 = zeros(nx,1);
pi_0(1) = 1;
pi_mc = zeros(nx,cf);
pi_mc(1,1) = 1;
for i = 1:cf-1
    pi_mc(:,i+1) = Pf*pi_mc(:,i);
end
pi_avg = (N0:-1:0)*pi_mc;

% Result Plot
figure(1); hold on % Time series graph
% plot(time_sim, Non)
plot(time_sim, mean(Non,2),'k-','linewidth', 2)
plot(t_mc, pi_avg,'r--')
xlabel('Day'); ylabel('# of sat.','linewidth', 2)

figure(2); hold on % Distribution at the last day
histogram(Non(end,:),'Normalization','probability')
plot((N0:-1:0), pi_mc(:,end)', 'r*')

