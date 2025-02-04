% This script vaildates the analysis method of direct spare management strategy.
% Proposed Markov Chain analysis results are compared with the Monte Carlo simulation results.

close all
clear all
clc

%.. PATH
mfilepath = pwd;
idcs = strfind(mfilepath,'\');
libdir = mfilepath(1:idcs(end));
addpath([libdir, 'CommonSource'])

%% Test Param
%.. Sim time
dt_sim = 1; % [day]
time_sim = 0:dt_sim:365*5000;

%.. Marcov Chain Period
dt_mc = dt_sim;  % [day]

%.. Failure rate
p_fail = 0.25/365; % [#/day] 0.05, 0.1, 0.15
p_sim = p_fail * dt_sim; % [#/dt_sim]
p_mc = p_fail * dt_mc; % [#/dt_mc]

%.. Failure type
p_type = 1; % 0:Constant failure, 1:State Dependent failure
n_sat = 40; % Number of operating satellite per orbit for nominal condition

%.. Direct (Q,R) Policy Parameter
Q = 4;
R = 40;

%.. State Parameter p.7
Xmax = Q + R; % Max State Level: bar(N_sat) 
Xnum = 0:1:Xmax; % State Counts: 0,1,...,bar(N_sat)

%.. Direct LV Parameters
mu_LV   =   40; % Mean lead time [day]
T0_LV   =   20; % Constant shift lead time [day]
dTlv_max = ceil((T0_LV + mu_LV + 5*mu_LV)/dt_mc); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]

%% Direct (Q,R) Policy 
rng('default')
iter_max = 1; % Number of different initial condition

% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock right after the resupply moment
Xq = zeros(Xmax+1, iter_max);

% The histogram of the  number of stock at the reordering moment
Xr =  zeros(Xmax+1, iter_max);

% Lead time distribution
Xlv = zeros(dTlv_max, iter_max); 

%% Run Each Simulation
for i = 1:iter_max
    % The number of satellite at current time step
    Non_k = R+round(Q*rand);        
    
    % ETC 
    cnt_lv = -1; % -1 for not ordered
    lv_cnt = 0; % # of LV order
    
    % Apply Policy
    for k = 1:length(time_sim)
        %%% 1. Generate Fail Sample at Every Contact
        if p_type == 0 %.. Const Failure Rate
            N_fail = CustomPoisRnd(n_sat*p_sim, 1);
        else %.. State Dependant Failure Rate
            if Non_k > n_sat
                N_fail = CustomPoisRnd(n_sat*p_sim, 1);
            else
                N_fail = CustomPoisRnd(Non_k*p_sim, 1);
            end
        end
        
        % Update the number of available stock
        Non_k = max([Non_k - N_fail, 0]);

        %%% 2. Check Resupply Arrival
        if cnt_lv == 0 % Arrive at this step
            % Update Non and Xq
            Non_k = Non_k + Q;
            Xq(Non_k+1,i) = Xq(Non_k+1,i) + 1; % +1 for index
            
            % Update LV Parameters
            cnt_lv = -1;
            lv_cnt = lv_cnt + 1;
        elseif cnt_lv > 0 % Wait for arrival
            cnt_lv = cnt_lv - 1;
        end

        %%% 3. Check Reorder at Every Time Step
        if cnt_lv == -1 && Non_k <= R % Reorder if not ordered and N <= R
            % Sample Lead Time and Save
            dT_LV = T0_LV + CustomExpRnd(mu_LV,1);
            dT_LV = ceil(dT_LV/dt_sim);
            if dTlv_max < dT_LV
                Xlv(end,i) = Xlv(end,i) + 1;
            else
                Xlv(dT_LV,i) = Xlv(dT_LV,i) + 1;
            end
            
            % Update Xr
            Xr(Non_k+1,i) = Xr(Non_k+1,i) + 1; % +1 for index
            
            % Save Remaining time step before arrival
            cnt_lv = dT_LV - 1;
        end
        
        % Save Stock Profile at current step after replinishment
        Non(k,i) = Non_k;
        Xon(Non_k+1,i) = Xon(Non_k+1,i) + 1; % +1 for index
    end
end

%% Plot Simulation Result
xx_edge = -0.5:1:(Xmax+0.5);
dT_edge = 0.5:1:(dTlv_max+0.5);

figure(1)
plot(time_sim, Non)

% figure(2); hold on
% histogram(Non(:),'Normalization','probability')
% xlabel('Number of stock for entire period')
% ylabel('Probability')

figure(2);
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xon,2),'Normalization','probability')
xlabel('Number of stock for entire period')
ylabel('Probability')

figure(3)
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xq,2),'Normalization','probability')
xlabel('Number of stock right after resupply')
ylabel('Probability')

figure(4)
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xr,2),'Normalization','probability')
xlabel('Number of stock at reordering of resupply')
ylabel('Probability')

figure(5); hold on
histogram(Non(Non > R),'Normalization','probability')
xlabel('Number of stock for non-reordering period')
ylabel('Probability')

figure(6); hold on
histogram(Non(Non <= R),'Normalization','probability')
xlabel('Number of stock for wait period')
ylabel('Probability')

figure(7)
histogram('BinEdges', dT_edge(1:end-1), 'BinCounts', sum(Xlv(1:end-1,:),2 ),'Normalization','probability')
xlabel('Lead Time Distribution')
ylabel('Probability')


%% Run Analysis Method
% As dt -> 0, pi_r becomes [0 0... 0 1 0... 0] 1 at N=R
f = p_mc;
lam = 1/mu_LV;
[xx, PI, T] = ExactDirectProb(p_fail,p_type,Q,R,lam,dt_mc,T0_LV,n_sat);


figure(2); hold on
plot(xx(1:13), PI.pi_dr(1:13), 'r*')
xlim([xx(13)-0.5, xx(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(3); hold on
plot(xx(1:13), PI.pi_q(1:13), 'r*')
xlim([xx(13)-0.5, xx(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(4); hold on
plot(xx(1:13), PI.pi_r(1:13), 'r*')
xlim([xx(13)-0.5, xx(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(5); hold on
plot(xx(1:13), PI.pi_np(1:13), 'r*')
xlim([xx(13)-0.5, xx(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(6); hold on
plot(xx(1:13), PI.pi_wp(1:13), 'r*')
xlim([xx(13)-0.5, xx(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

dtxx = T0_LV+1:dTlv_max;
figure(7); hold on
plot(dtxx, exp(-(dtxx - T0_LV - 1)/mu_LV)*(1-exp(-dt_mc/mu_LV)), 'r*')

%.. Given state distribution Compute Mean, ...
n_avg = sum(xx.*PI.pi_dr);
idx = find(xx == n_sat);
rho_loss = sum(PI.pi_dr(idx+1:end));
n_spare = sum((xx(1:idx-1)-n_sat).*PI.pi_dr(1:idx-1));

