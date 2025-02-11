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
iter_max = 10; % Number of different initial condition

%.. Sim time
dt_sim      =   1;                % [day]
time_sim    =   0:dt_sim:365*5000;

%.. Marcov Chain Period
dt_mc       =   dt_sim;             % [day]

%.. Number of Asymmetric Parking Orbit
n_park      =   4;

%.. In-plane Period (Time duration of in-plane for subsequent parking contact)
dt_plane    =   5*[13 15 25 40];     % [day] (Need n_park number of periods)
cnt_plane   =   round(dt_plane/dt_sim);
dT_plane    =   sum(dt_plane);      % Repeat period
Cnt_plane   =   sum(cnt_plane);     % Repeat period count

%.. Pre-compute the time when parking orbits ary synced.
Cnt_max = ceil(length(time_sim)/Cnt_plane);
rpt_cnt = repmat(cnt_plane',Cnt_max,1);
rpt_idx = repmat((1:n_park)',Cnt_max,1);
idx_contact = [0; cumsum(rpt_cnt)] + 1;
idx_contact(idx_contact > length(time_sim)) = [];
rpt_idx = rpt_idx(1:length(idx_contact));
cnt_park = zeros(length(time_sim),1);
cnt_park(idx_contact) = rpt_idx;

%.. Failure rate (test parameter, the demand distribution is not Poisson)
p_fail      =   0.5/365; % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_mc        =   p_fail * dt_mc;     % [#/dt_mc]

%.. Failure type
p_type = 0; % 0:Constant failure, 1:State Dependent failure
n_sat = 40; % Number of operating satellite per orbit for nominal condition

%.. Indirect In-plane (Q,R) Policy Parameter
Q       =   4;  % Assume same Q value for all parking orbit (shared spare-batch config.)
R       =   [37, 40, 43, 44]; % Apply different R value for each parking orbit 
%.. Rmk: To handle different (drift)lead time, allowing different R value will give better result.
%.. However, it will increase the decision variable a lot i.e., N_Park x N_Plane

%.. State Parameter
Xmax = Q + max(R); % Max State Level
Xnum = 0:1:Xmax; % State Counts

%.. Parking Availablity (Test distribution, each orbit will have different distribution)
Dmax = ceil(Xmax/Q); % Max Demend Level
Pav = sqrt(Dmax+1:-1:1); % Test parking stock distribution, Pav(k): Pr. having k-1 batch of spares
Pav = repelem(Pav,n_park,1) + 1.0*rand(n_park, length(Pav));
for i = 1:n_park
    Pav(i,:) = Pav(i,:)/sum(Pav(i,:));
end
Pav_sum = cumsum(Pav,2);
Kappa = [ones(n_park,1), 1 - Pav_sum(:,1:end-1)]; % Kappa(k): Pr. having more than k-1 batch of spares
Nav = 0:1:Dmax; % Demand Count

%% (Q, R) Policy with fixed time order
% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock right after the resupply moment
Xq = zeros(Xmax+1, n_park, iter_max);

% The histogram of the  number of stock at the reordering moment
Xr =  zeros(Xmax+1, n_park, iter_max);

% Parking available histogram
Xav = zeros(Dmax+1, n_park, iter_max);  
Xdmd = zeros(Dmax+1, n_park, iter_max); 

%% Run Each Simulation
for iter = 1:iter_max
    % The number of satellite at current time step
    Non_k = R+round(Q*rand);     

    % Apply Policy
    for i = 1:length(time_sim)
        % 1. Generate Fail Sample at Every Time Step
        if p_type == 0 %.. Const Failure Rate
            N_fail = CustomPoisRnd(n_sat*p_sim, 1);
        else %.. State Dependant Failure Rate
            if Non_k > n_sat
                N_fail = CustomPoisRnd(n_sat*p_sim, 1);
            else
                N_fail = CustomPoisRnd(Non_k*p_sim, 1);
            end
        end
        % Satellites after failure
        Non_k = max([Non_k - N_fail, 0]);

        % 2. Check Reorder at Every Lead Period
        if cnt_park(i)~= 0 % Contacted with the one of parking orbits
            idx_park = cnt_park(i); % Index of the contacted parking orbit 
            % Demand
            n_Req = R(idx_park) + 1 - Non_k;
            if n_Req > 0 % Reorder is required
                n_dmd = ceil(n_Req/Q);
                Xdmd(n_dmd+1,idx_park,iter) = Xdmd(n_dmd+1,idx_park,iter) + 1; % +1 for index
            else % Reorder is not required
                n_dmd = 0;
                Xdmd(1,idx_park,iter) = Xdmd(1,idx_park,iter) + 1;
            end
            
            % Update Xr
            Xr(Non_k+1,idx_park,iter) = Xr(Non_k+1,idx_park,iter) + 1; % +1 for index
            
            % Sample the number of available parking spares
            n_av = find(Pav_sum(idx_park,:) - rand >= 0, 1) - 1;
            Xav(n_av+1,idx_park,iter) = Xav(n_av+1,idx_park,iter) + 1; % +1 for index

            % Apply maximum feasible reorder #
            n_feas = min(n_av, n_dmd);
            Non_k = Non_k + n_feas*Q;
            Xq(Non_k+1,idx_park,iter) = Xq(Non_k+1,idx_park,iter) + 1; % +1 for index
        end

        % Save Stock Profile at current step after replinishment
        Non(i,iter) = Non_k;
        Xon(Non_k+1,iter) = Xon(Non_k+1,iter) + 1; % +1 for index
    end
end

%% Plot Simulation Result
close all

xx_edge = -0.5:1:(Xmax+0.5);
dmd_edge = -0.5:1:(Dmax+0.5);

figure(1)
plot(time_sim, Non)

figure(2);
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xon,2),'Normalization','probability')
xlabel('Number of stock for entire period')
ylabel('Probability')

figure(3)
histogram('BinEdges', xx_edge, 'BinCounts', sum(sum(Xq,3),2),'Normalization','probability')
xlabel('Number of stock right after resupply')
ylabel('Probability')

figure(31); hold on
tmp = sum(Xq,3);
for i = 1:n_park
    plot(Xnum,tmp(:,i)/sum(tmp(:,i)),'-o')
end
xlabel('Number of stock right after resupply')
ylabel('Probability')


figure(4)
histogram('BinEdges', xx_edge, 'BinCounts', sum(sum(Xr,3),2),'Normalization','probability')
xlabel('Number of stock at reordering of resupply')
ylabel('Probability')

figure(41); hold on
tmp = sum(Xr,3);
for i = 1:n_park
    plot(Xnum,tmp(:,i)/sum(tmp(:,i)),'-o')
end
xlabel('Number of stock at reordering of resupply')
ylabel('Probability')

figure(5); hold on
histogram('BinEdges', dmd_edge, 'BinCounts', sum(sum(Xdmd,3),2),'Normalization','probability')
xlabel('Number of demand')
ylabel('Probability')

figure(51); hold on
tmp = sum(Xdmd,3);
for i = 1:n_park
    plot(Nav,tmp(:,i)/sum(tmp(:,i)),'-o')
end
xlabel('Number of demand')
ylabel('Probability')

%% Run Analysis Method
[xx, PI, T] = AsymInDirectPlane(p_mc, p_type, n_sat, n_park, Kappa, Q, R, dt_mc, dt_plane);

figure(2); hold on
plot(xx, PI.pi_ir, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(3); hold on
plot(xx, mean(PI.pi_q,2), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(31); set(gca,'ColorOrderIndex',1)
for i = 1:n_park
    plot(xx,PI.pi_q(:,i),'*')
end

figure(4); hold on
plot(xx, mean(PI.pi_r,2), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(41); set(gca,'ColorOrderIndex',1)
for i = 1:n_park
    plot(xx,PI.pi_r(:,i),'*')
end

figure(5); hold on
plot(Nav, mean(PI.pi_dmd,2), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(51); set(gca,'ColorOrderIndex',1)
for i = 1:n_park
    plot(Nav, PI.pi_dmd(:,i), '*')
end



