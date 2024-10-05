% This script vaildates the analysis method of parking orbit under indirect spare management strategy.
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
dt_sim      =   0.5;                  % [day]
time_sim    =   0:dt_sim:365*1000;

%.. Marcov Chain Period
dt_mc       =   dt_sim;                 % [day]

%.. In-plane Period (Time duration of in-plane for subsequent parking contact, RAAN Drift time)
dt_plane    =   70;     % [day]
cnt_plane   =   round(dt_plane/dt_sim);

%.. Parking Period (Time duration of parking for subsequent in-plane contact)
dt_park  =   12;     % [day]
cnt_park =   round(dt_park/dt_mc);

%.. LV Lead Time (Modeled as bias shifted exponential distribution)
dt_mu_lv    =   50;     % [day]
dt_bias_lv  =   20;     % [day]
cnt_lv_max  =   ceil((dt_bias_lv + dt_mu_lv + 5*dt_mu_lv)/dt_mc); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]

%.. Parking (Q,R) Policy Parameter (For this, Q > R)
Qp      =   6;
Rp      =   3;

%.. State Parameter
Xmax = Qp + Rp; % Max State Level
Xnum = 0:1:Xmax; % State Counts

%.. Failure rate (test demand distribution)
p_fail      =   40*0.2/365; % [#/day]
p_drift     =   p_fail * dt_mc * cnt_plane;  % [#/dt_plane]


%% Parking (Q,R) Policy
rng('default')
iter_max = 5; % Number of different initial condition

% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock right after the resupply moment
Xq = zeros(Xmax+1, iter_max);

% The histogram of the  number of stock at the reordering moment
Xr =  zeros(Xmax+1, iter_max);

% Lead time distribution
Xnp = zeros(cnt_lv_max, iter_max); 
Xnp1 = zeros(cnt_lv_max, iter_max); 
Xlv = zeros(cnt_lv_max, iter_max); 

%% Run Each Simulation
for i = 1:iter_max
    % The number of satellite at current time step
    Non_k = Rp+round(Qp*rand);        
    
    % ETC 
    cnt_lv = -1; % -1 for not ordered
    cnt_np = 0;
    lv_cnt = 0;
    
    % Apply Policy
    for k = 1:length(time_sim)
        %%% 1. Check Resupply Arrival
        if cnt_lv == 0 % Arrive at this step
            % Update Non and Xq
            Non_k = Non_k + Qp;
            Xq(Non_k+1,i) = Xq(Non_k+1,i) + 1; % +1 for index
            
            % Update LV Parameters
            cnt_lv = -1;
            lv_cnt = lv_cnt + 1;

            % Reset Non-Reordering period count
            cnt_np = 0;
            cnt_np1 = ceil(k/cnt_park)*cnt_park - k;
            Xnp1(cnt_np1+1,i) = Xnp1(cnt_np1+1,i) + 1;
            

        elseif cnt_lv > 0 % Wait for arrival
            cnt_lv = cnt_lv - 1;

        elseif cnt_lv == -1 % Has not been ordered
            % Increase Non-Reordering period count
            cnt_np = cnt_np + 1;
        end

        %%% 2. Generate Fail Sample at Every Contact
        if mod(k,cnt_park) == 0
            N_dmd = poissrnd(p_drift, 1); %.. Const Failure Rate
        else
            N_dmd = 0;
        end
        
        % Update the number of available spares
        Non_k = max([Non_k - N_dmd, 0]);
        
        %%% 3. Check Reorder at Every Time Step
        if cnt_lv == -1 && Non_k <= Rp
            % Sample Lead Time and Save
            dt_lv = dt_bias_lv + exprnd(dt_mu_lv,1);
            cnt_lv = ceil(dt_lv/dt_mc); % If dt_lv: [0,1] -> 1

            if cnt_lv_max < cnt_lv
                Xlv(end,i) = Xlv(end,i) + 1;
            else
                Xlv(cnt_lv,i) = Xlv(cnt_lv,i) + 1;
            end
            
            % Update Xr
            Xr(Non_k+1,i) = Xr(Non_k+1,i) + 1; % +1 for index
            
            % Save Remaining time step before arrival
            cnt_lv = cnt_lv - 1;
        end
        
        % Save Stock Profile at current step after replinishment
        Non(k,i) = Non_k;
        Xon(Non_k+1,i) = Xon(Non_k+1,i) + 1; % +1 for index
    end
end

%% Plot Simulation Result
xx_edge = -0.5:1:(Xmax+0.5);
dt_edge = 0.5:1:(cnt_lv_max+0.5);

figure(1)
plot(time_sim, Non)

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
histogram(Non(Non > Rp),'Normalization','probability')
xlabel('Number of stock for non-reordering period')
ylabel('Probability')

figure(6); hold on
histogram(Non(Non <= Rp),'Normalization','probability')
xlabel('Number of stock for wait period')
ylabel('Probability')

figure(7)
histogram('BinEdges', dt_edge(1:end-1), 'BinCounts', sum(Xlv(1:end-1,:),2 ),'Normalization','probability')
xlabel('Lead Time Distribution')
ylabel('Probability')

figure(8)
histogram('BinEdges', dt_edge(1:end-1), 'BinCounts', sum(Xnp(1:end-1,:),2 ),'Normalization','probability')
xlabel('Lead Time Distribution')
ylabel('Probability')

%% Run Analysis Method
%.. In-plane demand parameter (test)
Eta = poisspdf(0:Xmax, p_drift);
[xx, PI, T] = ExactInDirectPark(Eta,Qp,Rp,dt_mc,dt_park,dt_mu_lv,dt_bias_lv);
pi_q = flip(sum(Xq,2)/sum(sum(Xq,2)));
pi_r = flip(sum(Xr,2)/sum(sum(Xr,2)));

figure(2); hold on
plot(xx, PI.pi_ir, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(3); hold on
plot(xx, PI.pi_q, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(4); hold on
plot(xx, PI.pi_r, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(5); hold on
plot(xx, PI.pi_np, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(6); hold on
plot(xx, PI.pi_wp, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

dtxx = dt_bias_lv:cnt_lv_max*dt_mc;
figure(7); hold on
plot(dtxx/dt_mc, exp(-(dtxx - dt_bias_lv)/dt_mu_lv)*(1-exp(-dt_mc/dt_mu_lv)),'r*')
