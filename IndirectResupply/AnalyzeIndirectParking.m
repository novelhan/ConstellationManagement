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
dt_sim      =   1;                  % [day]
time_sim    =   0:dt_sim:365*5000;

%.. In-plane Period (Time duration of in-plane for subsequent parking contact, RAAN Drift time)
dt_plane    =   70;     % [day]
cnt_plane   =   round(dt_plane/dt_sim);

%.. Parking Period (Time duration of parking for subsequent in-plane contact)
dt_park  =   15;     % [day]
cnt_park =   round(dt_park/dt_sim);

%.. LV Lead Time (Modeled as bias shifted exponential distribution)
mu_lv = 60;     % [day]
dt_lv = 20;     % [day]
cnt_lv_max = ceil((dt_lv + mu_lv + 5*mu_lv)/dt_sim); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]

%.. Parking (Q,R) Policy Parameter (For this, Q > R)
%.. Time based should only work for Q > R but it gives accurate result when Q = R. Need investigation.
Qp      =   8;
Rp      =   5;

%.. State Parameter
Xmax = Qp + Rp; % Max State Level
Xnum = 0:1:Xmax; % State Counts

%.. Failure rate (test demand distribution)
p_fail      =   40*0.2/365; % [#/day]
p_drift     =   p_fail * dt_plane;  % [#/dt_plane]


%% Parking (Q,R) Policy
rng('default')
iter_max = 5; % Number of different initial condition

% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock right after the resupply moment
Xq = zeros(Xmax+1, iter_max);

% The histogram of the number of stock at the reordering moment
Xr =  zeros(Xmax+1, iter_max);

% The histogram of the number of stock at the contact moment
Xc =  zeros(Xmax+1, iter_max);

% Lead time distribution
Xlv = zeros(cnt_lv_max, iter_max); 

%% Run Each Simulation
for i = 1:iter_max
    % The number of satellite at current time step
    Non_k = Rp+round(Qp*rand);        
    
    % ETC 
    cnt_lv = -1; % -1 for not ordered
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

        elseif cnt_lv > 0 % Wait for arrival
            cnt_lv = cnt_lv - 1;

        end

        %%% 2. Generate Fail Sample at Every Contact
        if mod(k,cnt_park) == 0
            N_dmd = CustomPoisRnd(p_drift, 1); %.. Const Failure Rate
            
            % Update Xc
            Xc(Non_k+1,i) = Xc(Non_k+1,i) + 1; % +1 for index
            
            % Update the number of available spares
            Non_k = max([Non_k - N_dmd, 0]);
        end
        
        %%% 3. Check Reorder at Every Time Step
        if cnt_lv == -1 && Non_k <= Rp
            % Sample Lead Time and Save
            dt_LV = dt_lv + CustomExpRnd(mu_lv,1);
            cnt_lv = ceil(dt_LV/dt_sim); % If dt_lv: [0,1] -> 1

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

figure(7); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xc,2),'Normalization','probability')
xlabel('Number of stock at the contact moment')
ylabel('Probability')

figure(8)
histogram('BinEdges', dt_edge(1:end-1), 'BinCounts', sum(Xlv(1:end-1,:),2 ),'Normalization','probability')
xlabel('Lead Time Distribution')
ylabel('Probability')

figure(9)
histogram('BinEdges', xx_edge, 'BinCounts', cumsum(sum(Xc,2),'reverse')/sum(sum(Xc,2)))
xlabel('Parking Spares Availability')
ylabel('Probability')


%% Run Analysis Method
%.. Marcov Chain Period
dt_mc       =   dt_sim;                 % [day]

%.. In-plane demand parameter (test)
Eta = CustomPoisPdf(0:Xmax, p_drift);

ParaParking.Q = Qp;
ParaParking.R = Rp;
ParaParking.dt_mc = dt_sim;
ParaParking.dt_park = dt_park;
ParaParking.mu_lv = mu_lv;
ParaParking.dt_lv = dt_lv;

%.. Method 1
ParaParking.method = 0;
disp('Time Based Method:')
tic
for i = 1:100
    [PI1, T1] = SolveInDirectPark(Eta, ParaParking);
end
toc

figure(2); hold on
plot(PI1.X, PI1.pi_avg, 'r*')

figure(3); hold on
plot(PI1.X, PI1.pi_q, 'r*')

figure(4); hold on
plot(PI1.X, PI1.pi_r, 'r*')

figure(5); hold on
plot(PI1.X, PI1.pi_np, 'r*')

figure(6); hold on
plot(PI1.X, PI1.pi_wp, 'r*')

figure(9); hold on
plot(Xnum, PI1.Pdav, 'r*')


%.. Method 2
ParaParking.method = 1;
disp('Ratio Based Method:')
tic
for i = 1:100
    [PI2, T2] = SolveInDirectPark(Eta, ParaParking);
end
toc

figure(2); hold on
plot(PI2.X, PI2.pi_avg, 'go')
legend('Sim.', 'Sol.1', 'Sol.2', 'location', 'best')

figure(3); hold on
plot(PI2.X, PI2.pi_q, 'go')
legend('Sim.', 'Sol.1', 'Sol.2', 'location', 'best')

figure(4); hold on
plot(PI2.X, PI2.pi_r, 'go')
legend('Sim.', 'Sol.1', 'Sol.2', 'location', 'best')

figure(5); hold on
plot(PI2.X, PI2.pi_np, 'go')
legend('Sim.', 'Sol.1', 'Sol.2', 'location', 'best')

figure(6); hold on
plot(PI2.X, PI2.pi_wp, 'go')
legend('Sim.', 'Sol.1', 'Sol.2', 'location', 'best')

figure(9); hold on
plot(Xnum, PI2.Pdav, 'go')
legend('Sim.', 'Sol.1', 'Sol.2', 'location', 'best')

dtxx = dt_lv:cnt_lv_max*dt_mc;
figure(8); hold on
plot(dtxx/dt_mc, exp(-(dtxx - dt_lv)/mu_lv)*(1-exp(-dt_mc/mu_lv)),'r*')

disp('Simulation based avg. cycle period')
disp(time_sim(end)*iter_max/sum(sum(Xlv,1)))
disp('')

disp('Analysis based avg. cycle period')
disp([T1.T_avg, T2.T_avg])
disp('')

