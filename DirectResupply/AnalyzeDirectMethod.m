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
rng('default')
iter_max = 5; % Number of different initial condition

%.. Sim time
dt_sim = 1; % [day]
time_sim = 0:dt_sim:365*100;

%.. Failure rate
p_fail = 0.1/365; % [#/day] 0.05, 0.1, 0.15
p_sim = p_fail * dt_sim; % [#/dt_sim]

%.. Failure type
p_type = 1; % 0:Constant failure, 1:State Dependent failure
N_sat = 40; % Number of operating satellite per orbit for nominal condition

%.. Direct (Q,R) Policy Parameter
Q = 3;
R = 40;

%.. State Parameter p.7
Xmax = Q + R; % Max State Level: bar(N_sat) 
Xnum = 0:1:Xmax; % State Counts: 0,1,...,bar(N_sat)

%.. Direct LV Parameters
mu_LV   =   40; % Mean lead time [day]
T0_LV   =   20; % Constant shift lead time [day]
dTlv_max = ceil((T0_LV + mu_LV + 5*mu_LV)/dt_sim); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]

%% Direct (Q,R) Policy 
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
for itr = 1:iter_max
    % The number of satellite at current time step
    Non_k = R+round(Q*rand);        
    
    % ETC 
    cnt_lv = -1; % -1 for not ordered
    lv_cnt = 0; % # of LV order
    
    % Apply Policy
    for k = 1:length(time_sim)
        %%% 1. Generate Fail Sample at Every Contact
        if p_type == 0 %.. Const Failure Rate
            N_fail = CustomPoisRnd(N_sat*p_sim, 1);
        else %.. State Dependant Failure Rate
            if Non_k > N_sat
                N_fail = CustomPoisRnd(N_sat*p_sim, 1);
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
            Xq(Non_k+1,itr) = Xq(Non_k+1,itr) + 1; % +1 for index
            
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
                Xlv(end,itr) = Xlv(end,itr) + 1;
            else
                Xlv(dT_LV,itr) = Xlv(dT_LV,itr) + 1;
            end
            
            % Update Xr
            Xr(Non_k+1,itr) = Xr(Non_k+1,itr) + 1; % +1 for index
            
            % Save Remaining time step before arrival
            cnt_lv = dT_LV - 1;
        end
        
        % Save Stock Profile at current step after replinishment
        Non(k,itr) = Non_k;
        Xon(Non_k+1,itr) = Xon(Non_k+1,itr) + 1; % +1 for index
    end
end

%% Plot Simulation Result
xx_edge = -0.5:1:(Xmax+0.5);
dT_edge = 0.5:1:(dTlv_max+0.5);
Xon = sum(Xon,2);
Xq = sum(Xq,2);
Xr = sum(Xr,2);

figure(1)
plot(time_sim, Non)

% figure(2); hold on
% histogram(Non(:),'Normalization','probability')
% xlabel('Number of stock for entire period')
% ylabel('Probability')

figure(2);
histogram('BinEdges', xx_edge, 'BinCounts', Xon,'Normalization','probability')
xlabel('Number of stock for entire period')
ylabel('Probability')

figure(3)
histogram('BinEdges', xx_edge, 'BinCounts', Xq,'Normalization','probability')
xlabel('Number of stock right after resupply')
ylabel('Probability')

figure(4)
histogram('BinEdges', xx_edge, 'BinCounts', Xr,'Normalization','probability')
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

% figure(7)
% histogram('BinEdges', dT_edge(1:end-1), 'BinCounts', sum(Xlv(1:end-1,:),2 ),'Normalization','probability')
% xlabel('Lead Time Distribution')
% ylabel('Probability')

%% Run Analysis Method
ParaInPlane.f_ref = p_fail;
ParaInPlane.f_type = p_type;
ParaInPlane.N_sat = N_sat;
ParaInPlane.Q = Q;
ParaInPlane.R = R;
ParaInPlane.dt_mc = dt_sim/2;
ParaInPlane.mu_lv = mu_LV;
ParaInPlane.dt_lv = T0_LV;

%.. Method 1
disp('Time Based Method:')
ParaInPlane.method = 0;
tic
for i = 1:100
[PI1, T1] = SolveDirectProb(ParaInPlane);
end
toc 

figure(2); hold on
plot(PI1.X, PI1.pi_dr, 'r*')

figure(3); hold on
plot(PI1.X, PI1.pi_q, 'r*')

figure(4); hold on
plot(PI1.X, PI1.pi_r, 'r*')

figure(5); hold on
plot(PI1.X, PI1.pi_np, 'r*')

figure(6); hold on
plot(PI1.X, PI1.pi_wp, 'r*')

%.. Method 2
disp('Ratio Based Method:')
ParaInPlane.method = 1;
tic
for i = 1:100
    [PI2, T2] = SolveDirectProb(ParaInPlane);
end
toc

figure(2); hold on
plot(PI2.X, PI2.pi_dr, 'go')
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

% dtxx = T0_LV+1:dTlv_max;
% figure(7); hold on
% plot(dtxx, exp(-(dtxx - T0_LV - 1)/mu_LV)*(1-exp(-dt_mc/mu_LV)), 'r*')

%% Additional Cost Validation
%.. Cost Model
c_build = 0.5; % Manufacturing cost per satellite [M$/sat]
c_hold = 0.5/365; % Holding cost per satellite in In-Plane orbit [M$/sat/day]
c_lv_small_part = 3; %[M$/sat/launch] for direct LV (small size)
c_lv_small_full = 7.5;%[M$/launch] for direct LV (small size)
Qmax = 3;

%.. Simulation Result
Son = Non - N_sat;
Son(Son < 0) = 0; 
Navg = mean(mean(Non,1));
Savg = mean(mean(Son,1));
Nlv = sum(sum(Xlv,2));
f_lv = Nlv/(iter_max*time_sim(end));
f_sat = Q*Nlv/(iter_max*time_sim(end));
p_loss = sum(Xon(1:N_sat))/sum(Xon);

disp(' ')
disp('Simulation Results')
disp(['Avg. # Sat: ',num2str(Navg)])
disp(['Avg. # Spares: ',num2str(Savg)])
disp(['Build Cost: ',num2str(c_build*f_sat)])
disp(['Holding Cost: ',num2str(c_hold*Savg)])
if Q == Qmax
    disp(['Launch Cost: ',num2str(c_lv_small_full*f_lv)])
else
    disp(['Launch Cost: ',num2str(c_lv_small_part*Q*f_lv)])
end
disp(['P(Xi< N_sat): ',num2str(p_loss)])

%.. Analysis result
n_avg = sum(PI2.X.*PI2.pi_dr);
idx = find(PI2.X == N_sat);
p_loss = sum(PI2.pi_dr(idx+1:end));
n_spare = sum((PI2.X(1:idx-1)-N_sat).*PI2.pi_dr(1:idx-1));

disp(' ')
disp('Analysis Results')
disp(['Avg. # Sat: ', num2str(n_avg)])
disp(['Avg. # Spares: ',num2str(n_spare)])
disp(['Build Cost: ',num2str(c_build*Q/T2.T_dr)])
disp(['Holding Cost: ',num2str(c_hold*n_spare)])
if Q == Qmax
    disp(['Launch Cost: ',num2str(c_lv_small_full/T2.T_dr)])
else
    disp(['Launch Cost: ',num2str(c_lv_small_part*Q/T2.T_dr)])
end
disp(['P(Xi< N_sat): ',num2str(p_loss)])

%.. Using external function
ParaInPlane.N_plane = 1;
ParaCost.c_build = c_build;
ParaCost.c_hold = c_hold;
ParaCost.c_lv_small_part = c_lv_small_part;
ParaCost.c_lv_small_full = c_lv_small_full;
ParaConst.p_loss = 0;
[J, Cost] = CostDirectResupply([Q, R], ParaCost, ParaInPlane);
c = ConstDirectResupply([Q, R], ParaConst, ParaInPlane);

disp(' ')
disp('External Function Results')
disp(['Build Cost: ',num2str(Cost.C_build)])
disp(['Holding Cost: ',num2str(Cost.C_hold)])
disp(['Launch Cost: ',num2str(Cost.C_launch)])
disp(['P(Xi< N_sat): ',num2str(c)])
