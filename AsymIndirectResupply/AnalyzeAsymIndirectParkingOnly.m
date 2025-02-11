% This script vaildates the analysis method of parking orbit under indirect spare management
% strategy with assumption of asymmetric parking orbit distribution. Proposed Markov Chain analysis 
% results are compared with the Monte Carlo simulation results.

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
time_sim    =   0:dt_sim:365*5000;

%.. Marcov Chain Period
dt_mc       =   dt_sim;                 % [day]

%.. Number of Asymmetric Constellation Orbit
n_plane     =   4;

%.. Parking Period (Time duration of parking for subsequent in-plane contact)
%.. dt_park(i): Time to meet (i+1)th in-plane since the (i)th in-plane contact
dt_park  =   3*[3 6 9 12];     % [day] (Need n_plane number of periods)
cnt_park =   round(dt_park/dt_mc);
dT_park = sum(dt_park); % Repeat period (After this much time, the process repeat again)
Cnt_park = sum(cnt_park); % Repeat period count

%.. Pre-compute the time when in-plane orbits ary synced.
Cnt_max = ceil(length(time_sim)/Cnt_park);
rpt_cnt = repmat(cnt_park',Cnt_max,1);
rpt_idx = repmat((1:n_plane)',Cnt_max,1);
idx_contact = [0; cumsum(rpt_cnt)] + 1;
idx_contact(idx_contact > length(time_sim)) = [];
rpt_idx = rpt_idx(1:length(idx_contact));
cnt_plane = zeros(length(time_sim),1);
cnt_plane(idx_contact) = rpt_idx;

%.. LV Lead Time (Modeled as bias shifted exponential distribution)
dt_mu_lv    =   25;     % [day]
dt_bias_lv  =   10;     % [day]
cnt_lv_max  =   ceil((dt_bias_lv + dt_mu_lv + 5*dt_mu_lv)/dt_mc); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]

%.. Parking (Q,R) Policy Parameter (Q > R condition is removed by applying the second method)
Qp      =   2;
Rp      =   [1 2 3 4];
%.. Rmk: Unlike In-Plane transfer, single R value coulb be used. This is because, when to use 
%   different R value is not obvious. To have a general expression, we used multiple R values and
%   assumed different R value is used in between before and upcoming contact point.
%   Ex. We are in - 1 --- now --- 2 - : then R(1) is used.
%   But in the end, I think it is better to consider online optimization which is structurally 
%   unconstrained policy free approach. e.g., Determine when to order for given state.

%.. State Parameter
Xmax = Qp + max(Rp); % Max State Level
Xnum = 0:1:Xmax; % State Counts

%.. In-Plane Spare Demand (test parameter, simply used the Poisson with different input param.)
%.. Rmk: The demand distribution of each in-plane orbit is not Poisson and different each other
p_fail = 40*0.1/365; % [#/day]
p_drift = p_fail * dt_mc * dt_park;  % [#/dt_plane]
Eta = zeros(n_plane, Xmax+1);
for i = 1:n_plane
    Eta(i,:) = CustomPoisPdf(0:Xmax, p_drift(i));
end


%% Parking (Q,R) Policy
% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock right after the resupply moment
Xq = zeros(Xmax+1, iter_max);

% The histogram of the  number of stock at the reordering moment
Xr =  zeros(Xmax+1, n_plane, iter_max);

% The histogram of the number of stock at the contact moment
Xc =  zeros(Xmax+1, n_plane, iter_max);

%
Xp = zeros(length(time_sim), iter_max);

% Lead time distribution
Xlv = zeros(cnt_lv_max, iter_max); 

%% Run Each Simulation
for i = 1:iter_max
    % The number of satellite at current time step
    Non_k = Rp(1)+round(Qp*rand);        
    
    % ETC 
    cnt_lv = -1; % -1 for not ordered
    lv_cnt = 0;
    tmp = find(cnt_plane ~= 0,1);
    idx_plane = cnt_plane(tmp); % Find the index of previously contacted In-Plane 
    
    % Apply Policy
    for k = 1:length(time_sim)
        %% 1. Check Resupply Arrival
        if cnt_lv == 0 % Arrive at this step
            % Update Non and Xq
            Non_k = Non_k + Qp;
            Xq(Non_k+1,i) = Xq(Non_k+1,i) + 1; % +1 for index

            % Update LV Parameters
            cnt_lv = -1;

        elseif cnt_lv > 0 % Wait for arrival
            cnt_lv = cnt_lv - 1;

        elseif cnt_lv == -1 % Has not been ordered
            Xp(k,i) = 1;

        end

        %% 2. Send spares to the In-Plane at Every Contact moment
        if cnt_plane(k) ~= 0
            %.. Generate Demand of the contacted In-Plane
            idx_plane = cnt_plane(k);
            
            % Update Xc
            Xc(Non_k+1,idx_plane,i) = Xc(Non_k+1,idx_plane,i) + 1; % +1 for index
            
            % Eta(1):Pr. of 0 dmd, Eta(i):Pr. of i-1 dmd.
            N_dmd = find(rand < cumsum(Eta(idx_plane,:)), 1) - 1;
            
            % Update the number of available spares
            Non_k = max([Non_k - N_dmd, 0]); % Send every spares it has
        end
        
        
        %% 3. Check Reorder at Every Time Step
        if cnt_lv == -1 && Non_k <= Rp(idx_plane)
            % Sample Lead Time and Save
            dt_lv = dt_bias_lv + exprnd(dt_mu_lv,1);
            cnt_lv = ceil(dt_lv/dt_mc); % If dt_lv: [0,1] -> 1

            if cnt_lv_max < cnt_lv
                Xlv(end,i) = Xlv(end,i) + 1;
            else
                Xlv(cnt_lv,i) = Xlv(cnt_lv,i) + 1;
            end
            lv_cnt = lv_cnt + 1;
            
            % Update Xr
            Xr(Non_k+1,idx_plane,i) = Xr(Non_k+1,idx_plane,i) + 1; % +1 for index
            
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
histogram('BinEdges', xx_edge, 'BinCounts', sum(sum(Xr,3),2),'Normalization','probability')
xlabel('Number of stock at reordering of resupply')
ylabel('Probability')

figure(41); hold on
tmp = sum(Xr,3);
for i = 1:n_plane
    plot(Xnum,tmp(:,i)/sum(tmp(:,i)),'-o')
end
xlabel('Number of stock at reordering of resupply')
ylabel('Probability')

disp( 'Reorder Period Ratio:')
disp( sum(sum(Xr,3),1)/ sum(sum(sum(Xr,3),1)))

% figure(5); hold on
% histogram(Non(Non > Rp),'Normalization','probability')
% xlabel('Number of stock for non-reordering period')
% ylabel('Probability')
% 
% figure(6); hold on
% histogram(Non(Non <= Rp),'Normalization','probability')
% xlabel('Number of stock for wait period')
% ylabel('Probability')

figure(7)
histogram('BinEdges', dt_edge(1:end-1), 'BinCounts', sum(Xlv(1:end-1,:),2 ),'Normalization','probability')
xlabel('Lead Time Distribution')
ylabel('Probability')


%% Run Analysis Method
%.. In-plane demand parameter (test)
tic
[xx, PI, T] = AsymInDirectPark(Eta,Qp,Rp,n_plane,dt_mc,dt_park,dt_mu_lv,dt_bias_lv);
toc

figure(2); hold on
plot(xx, PI.pi_ir, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(3); hold on
plot(xx, sum(PI.pi_q,2)/sum(sum(PI.pi_q,2)), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(4); hold on
plot(xx, sum(PI.pi_r,2)/sum(sum(PI.pi_r,2)), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(41); set(gca,'ColorOrderIndex',1)
for i = 1:n_plane
    plot(xx,PI.pi_r(:,i)/sum(PI.pi_r(:,i)),'*')
end


%%%%
nx = length(xx);
[~,tt] = size(PI.PI_full);

figure(); hold on
plot(1:tt, PI.pi_full,'o')

Xtmp = zeros(nx,tt);
for i = 1:tt
    for n = 1:nx
        Xtmp(n,i) = sum(sum(Non(i:tt:length(Non),:) == (n-1)));
    end
end
Xtmp = Xtmp/sum(Xtmp(:,1));
set(gca,'ColorOrderIndex',1)
plot(1:tt, flip(Xtmp),'-')

%%%

ratio = zeros(2,tt);
for i = 1:tt
    ratio(1,i) = sum(PI.PI_full(1:nx,i));
    ratio(2,i) = sum(PI.PI_full(1+nx:end,i));
%     ratio(3,i) = 1 - ratio(1,i) - ratio(2,i);
end

figure(); hold on
plot(1:tt, ratio,'o')

Xtmp = zeros(2,tt);
for i = 1:tt
    for n = 1:2
        Xtmp(n,i) = sum(sum(Xp(i:tt:length(Non),:) == (n-1)));
    end
end
Xtmp = Xtmp/sum(Xtmp(:,1));

set(gca,'ColorOrderIndex',1)
plot(1:tt, flip(Xtmp),'-')



% 
% 
% dtxx = dt_bias_lv:cnt_lv_max*dt_mc;
% figure(7); hold on
% plot(dtxx/dt_mc, exp(-(dtxx - dt_bias_lv)/dt_mu_lv)*(1-exp(-dt_mc/dt_mu_lv)),'r*')

disp('Parking Availability')
disp(cumsum(sum(Xc,3),'reverse')./sum(sum(Xc,3)))
