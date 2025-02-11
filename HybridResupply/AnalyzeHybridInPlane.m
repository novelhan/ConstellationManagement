% This script vaildates the analysis method of in-plane orbit under the hybrid spare management strategy.
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
dt_sim      =   1;                  % [day]
time_sim    =   0:dt_sim:365*1000;

%.. Marcov Chain Period
dt_mc       =   1;             % [day]

%.. In-plane Period (Time duration of in-plane for subsequent parking contact, RAAN Drift time)
dt_plane   =   60;                 % [day]
cnt_plane  =   round(dt_plane/dt_sim);

%.. Failure rate (test parameter)
p_fail  =   0.4/365; % [#/day] 0.05, 0.1, 0.15
p_sim   =   p_fail * dt_sim;    % [#/dt_sim]

%.. Failure type
p_type = 0; % 0:Constant failure, 1:State Dependent failure
n_sat = 40; % Number of operating satellite per orbit for nominal condition

%.. Hybrid In-plane (Q1,R1,Q2,R2) Policy Parameter
Q1  =   4;
R1  =   40;
Q2  =   2;
R2  =   40;

%.. State Parameter
Xmax = Q1+Q2+max(R1,R2); % Max State Level
Xnum = 0:1:Xmax; % State Counts

%.. Parking Availablity (Test distribution)
Dmax = ceil(Xmax/Q1); % Max Demend Level
Pav = sqrt(Dmax+1:-1:1); % Test parking stock distribution, Pav(k): Pr. having k-1 batch of spares
Pav = Pav/sum(Pav); 
Pav_sum = cumsum(Pav);
Kappa = [1, 1 - Pav_sum(1:end-1)]; % Kappa(k): Pr. having more than k-1 batch of spares
Nav = 0:1:Dmax; % Demand Count

%.. Direct LV Parameters (Modeled as bias shifted exponential distribution)
mu_LV   =   30;     % [day]
dt_LV   =   15;
dTlv_max = ceil(dt_LV + mu_LV + 5*mu_LV); % Max Lead Time Bin (mean + 5*sigma) [day]
dTlv = 1:dTlv_max; % Lead Time Bin (Minimum is set as dt_sim)

%% (Q1, R1, Q2, R2) Policy
% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock for the indirect resupply moment
Xq1 = zeros(Xmax+1, iter_max);
Xr1 = zeros(Xmax+1, iter_max);

% The histogram of the number of stock for the direct resupply moment
Xq2 = zeros(Xmax+1, cnt_plane, iter_max);
Xr2 = zeros(Xmax+1, cnt_plane, iter_max);

% Direct Lead time distribution
Xlv = zeros(dTlv_max, iter_max); 

% Parking available histogram
Xav = zeros(Dmax+1, iter_max);  
Xdmd = zeros(Dmax+1, iter_max); 

%% Run Each Simulation
for itr = 1:iter_max
    % The number of satellite at current time step
    Non_k = Xmax - round(Q1*rand) - round(Q2*rand);
    
    % ETC 
    cnt_lv = -1; % -1 for not ordered
    lv_cnt = 0;
    
    % Apply Policy
    for k = 1:length(time_sim)
        % RAAN Contact Counter
        cnt_p = mod(k,cnt_plane) + 1; % 1, 2, ... , cnt_plane
        
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
        
        %%% 2. Check Direct Resupply Arrival
        if cnt_lv == 0 % Arrive at this step
            % Update Non and Xq
            Non_k = Non_k + Q2;
            Xq2(Non_k+1,cnt_p,itr) = Xq2(Non_k+1,cnt_p,itr) + 1; % +1 for index
            
            % Update LV Parameters
            cnt_lv = -1;
            lv_cnt = lv_cnt + 1;

        elseif cnt_lv > 0 % Wait for arrival
            cnt_lv = cnt_lv - 1;

        end
        
        %%% 3. Check Q1 Resupply when it contact with one of the parking orbit
        if cnt_p == 1
            % Demand
            n_Req = R1 + 1 - Non_k;
            if n_Req > 0 % Reorder is required
                n_dmd = ceil(n_Req/Q1);
                Xdmd(n_dmd+1,itr) = Xdmd(n_dmd+1,itr) + 1; % +1 for index
            else % Reorder is not required
                n_dmd = 0;
                Xdmd(1,itr) = Xdmd(1,itr) + 1;
            end
            
            % Update Xr
            Xr1(Non_k+1,itr) = Xr1(Non_k+1,itr) + 1; % +1 for index
            
            % Sample the number of available parking spares
            n_park = find(Pav_sum - rand >= 0, 1) - 1;
            Xav(n_park+1,itr) = Xav(n_park+1,itr) + 1; % +1 for index
            
            % Apply maximum feasible reorder #
            n_feas = min(n_park, n_dmd);
            Non_k = Non_k + n_feas*Q1;
            Xq1(Non_k+1,itr) = Xq1(Non_k+1,itr) + 1; % +1 for index
        end
        
        %%% 4. Check Q2 Resupply at Every Time Step
        if cnt_lv == -1 && Non_k <= R2
            % Sample Lead Time and Save
            dT_LV = dt_LV + CustomExpRnd(mu_LV,1);
            dT_LV = ceil(dT_LV/dt_sim);

            if dTlv_max < dT_LV
                Xlv(end,itr) = Xlv(end,itr) + 1;
            else
                Xlv(dT_LV,itr) = Xlv(dT_LV,itr) + 1;
            end
            
            % Update Xr
            Xr2(Non_k+1,cnt_p,itr) = Xr2(Non_k+1,cnt_p,itr) + 1; % +1 for index
            
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
dmd_edge = -0.5:1:(Dmax+0.5);

figure(1)
plot(time_sim, Non)

figure(2); hold on
histogram(Non(:),'Normalization','probability')
xlabel('Number of stock for entire period')
ylabel('Probability')

figure(3); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xq1,2), 'Normalization','probability')
xlabel('Number of stock right after Q1 resupply')
ylabel('Probability')

figure(4); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xr1,2), 'Normalization','probability')
xlabel('Number of stock at reordering of R1 resupply')
ylabel('Probability')

figure(5); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(sum(Xq2,3),2), 'Normalization','probability')
xlabel('Number of stock right after Q2 resupply')
ylabel('Probability')

figure(6); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(sum(Xr2,3),2), 'Normalization','probability')
xlabel('Number of stock at reordering of R2 resupply')
ylabel('Probability')

figure(7); hold on
histogram('BinEdges', dmd_edge, 'BinCounts', sum(Xdmd,2),'Normalization','probability')
xlabel('Number of demand')
ylabel('Probability')

Pbq2 = sum(Xq2,3);
Pbr2 = sum(Xr2,3);
for i = 1:cnt_plane
    Pbq2(:,i) = Pbq2(:,i)/sum(Pbq2(:,i));
    Pbr2(:,i) = Pbr2(:,i)/sum(Pbr2(:,i));
end

figure(8); hold on
for m = 1:(R2+1)
    plot(1:cnt_plane, [Pbr2(m,2:end), Pbr2(m,1)])
end
xlabel('Count after contact with parking orbit')
ylabel('Prob. of \pi_{r2} over time')

figure(9); hold on
for m = 1:(Xmax+1)
    plot(1:cnt_plane, [Pbq2(m,2:end), Pbq2(m,1)])
end
xlabel('Count after contact with parking orbit')
ylabel('Prob. of \pi_{q2} over time')

disp('Ratio between Q1 review period and Q2 order period is:')
disp(['Sim Result: ', num2str(mean(sum(Xlv)./sum(Xq1)))])

%% Run Analysis Code
ParaFail.dt_mc = dt_mc;
ParaFail.f_ref = p_fail;
ParaFail.f_type = p_type;
ParaFail.n_sat = n_sat;

%.. Direct Resupply Parameter
ParaDirect.mu = mu_LV;
ParaDirect.dt_lv = dt_LV;
ParaDirect.Q2 = Q2;
ParaDirect.R2 = R2;

%.. Indirect Resupply Parameter
ParaIndirect.dt_plane = dt_plane;
ParaIndirect.kappa = Kappa;
ParaIndirect.Q1 = Q1;
ParaIndirect.R1 = R1;

%.. Full state result
ParaDim.flag = 0;
tic
for i = 1:20
[x, PI, T] = HybridInPlane(ParaFail, ParaDirect, ParaIndirect, ParaDim);
end
toc
disp(['Full Result: ', num2str(T.T_q1/T.T_q2)])

s = 'r*';
figure(2); hold on
plot(x, mean(PI.pi_hr,2), s)

figure(3); hold on
plot(x, PI.pi_q1, s)

figure(4); hold on
plot(x, PI.pi_r1, s)

figure(5); hold on
plot(x, PI.pi_q2, s)

figure(6); hold on
plot(x, PI.pi_r2, s)

figure(7); hold on
plot(0:length(PI.pi_dmd)-1, PI.pi_dmd, s)

figure(8); hold on
set(gca,'ColorOrderIndex',1)
PI_r2 = flip(PI.PI_r2);
for m = 1:length(x)
    plot((1:round(dt_plane/dt_mc))*dt_mc, [PI_r2(m,2:end),PI_r2(m,1)], 'o')
end

PI_q2 = flip(PI.PI_q2);
figure(9); hold on
set(gca,'ColorOrderIndex',1)
for m = 1:length(x)
    plot((1:round(dt_plane/dt_mc))*dt_mc, [PI_q2(m,2:end), PI_q2(m,1)], 'o')
end

figure(10); hold on
set(gca,'ColorOrderIndex',1)
for m = 1:10
    plot((1:round(dt_plane/dt_mc))*dt_mc, PI.pi_hr(m,:), 'o')
end



%.. Reduced state result
ParaDim.flag = 1;
ParaDim.x_min = min(R1,R2) - ceil((2.5*p_sim*n_sat)*(dt_LV + 1.5*mu_LV)); % 2 sigma for both failure and lead time
ParaDim.n_seg = round( (dt_LV/dt_mc - 1)/3 );
tic
for i = 1:20
[x, PI, T] = HybridInPlane(ParaFail, ParaDirect, ParaIndirect, ParaDim);
end
toc
disp(['Reduced Result: ', num2str(T.T_q1/T.T_q2)])
s = 'mx';
figure(2); hold on
plot(x, mean(PI.pi_hr,2), s)
legend('Sim.', 'Full Sol.', 'Reduced Sol.', 'location', 'best')

figure(3); hold on
plot(x, PI.pi_q1, s)
legend('Sim.', 'Full Sol.', 'Reduced Sol.', 'location', 'best')

figure(4); hold on
plot(x, PI.pi_r1, s)
legend('Sim.', 'Full Sol.', 'Reduced Sol.', 'location', 'best')

figure(5); hold on
plot(x, PI.pi_q2, s)
legend('Sim.', 'Full Sol.', 'Reduced Sol.', 'location', 'best')

figure(6); hold on
plot(x, PI.pi_r2, s)
legend('Sim.', 'Full Sol.', 'Reduced Sol.', 'location', 'best')

figure(7); hold on
plot(0:length(PI.pi_dmd)-1, PI.pi_dmd, s)
legend('Sim.', 'Full Sol.', 'Reduced Sol.', 'location', 'best')

[~, b] = size(PI.PI_q2);

figure(8); hold on
set(gca,'ColorOrderIndex',1)
PI_r2 = flip([PI.PI_r2; zeros(x(end),b)]);
[a, ~] = size(PI_r2);
for m = 1:a
    plot((1:round(dt_plane/dt_mc))*dt_mc, [PI_r2(m,2:end),PI_r2(m,1)], 'x')
end

PI_q2 = flip([PI.PI_q2; zeros(x(end),b)]);
figure(9); hold on
set(gca,'ColorOrderIndex',1)
for m = 1:a
    plot((1:round(dt_plane/dt_mc))*dt_mc, [PI_q2(m,2:end), PI_q2(m,1)], 'x')
end

figure(10); hold on
set(gca,'ColorOrderIndex',1)
for m = 1:10
    plot((1:round(dt_plane/dt_mc))*dt_mc, PI.pi_hr(m,:), 'x')
end



