% This script vaildates the analysis method of in-plane orbit under the indirect spare management strategy.
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
dt_sim      =   1;                % [day]
time_sim    =   0:dt_sim:365*1000;

%.. Marcov Chain Period
dt_mc       =   dt_sim;             % [day]

%.. In-plane Period (Time duration of in-plane for subsequent parking contact, RAAN Drift time)
dt_plane    =   20;     % [day]
cnt_plane   =   round(dt_plane/dt_sim);

%.. Failure rate (test parameter)
p_fail  =   0.05/365; % [#/day] 0.05, 0.1, 0.15
p_sim   =   p_fail * dt_sim;    % [#/dt_sim]

%.. Failure type
p_type = 0; % 0:Constant failure, 1:State Dependent failure
n_sat = 40; % Number of operating satellite per orbit for nominal condition

%.. Indirect In-plane (Q,R) Policy Parameter
Q       =   4;
R       =   42;

%.. State Parameter
Xmax = Q + R; % Max State Level
Xnum = 0:1:Xmax; % State Counts

%.. Parking Availablity (Test distribution)
Dmax = ceil(Xmax/Q); % Max Demend Level
Pav = sqrt(Dmax+1:-1:1); % Test parking stock distribution, Pav(k): Pr. having k-1 batch of spares
Pav = Pav/sum(Pav); 
Pav_sum = cumsum(Pav);
Kappa = [1, 1 - Pav_sum(1:end-1)]; % Kappa(k): Pr. having more than k-1 batch of spares
Nav = 0:1:Dmax; % Demand Count

%% (Q, R) Policy with fixed time order
rng('default')
iter_max = 10; % Number of different initial condition

% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock right after the resupply moment
Xq = zeros(Xmax+1, iter_max);

% The histogram of the  number of stock at the reordering moment
Xr =  zeros(Xmax+1, iter_max);

% Parking available histogram
Xav = zeros(Dmax+1, iter_max);  
Xdmd = zeros(Dmax+1, iter_max); 

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
        if mod(i,cnt_plane) == 0
            % Demand
            n_Req = R + 1 - Non_k;
            if n_Req > 0 % Reorder is required
                n_dmd = ceil(n_Req/Q);
                Xdmd(n_dmd+1,iter) = Xdmd(n_dmd+1,iter) + 1; % +1 for index
            else % Reorder is not required
                n_dmd = 0;
                Xdmd(1,iter) = Xdmd(1,iter) + 1;
            end
            
            % Update Xr
            Xr(Non_k+1,iter) = Xr(Non_k+1,iter) + 1; % +1 for index
            
            % Sample the number of available parking spares
            n_park = find(Pav_sum - rand >= 0, 1) - 1;
            Xav(n_park+1,iter) = Xav(n_park+1,iter) + 1; % +1 for index

            % Apply maximum feasible reorder #
            n_feas = min(n_park, n_dmd);
            Non_k = Non_k + n_feas*Q;
            Xq(Non_k+1,iter) = Xq(Non_k+1,iter) + 1; % +1 for index

        end

        % Save Stock Profile at current step after replinishment
        Non(i,iter) = Non_k;
        Xon(Non_k+1,iter) = Xon(Non_k+1,iter) + 1; % +1 for index
    end
end

%% Plot Simulation Result
xx_edge = -0.5:1:(Xmax+0.5);
dmd_edge = -0.5:1:(Dmax+0.5);

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
histogram('BinEdges', dmd_edge, 'BinCounts', sum(Xdmd,2),'Normalization','probability')
xlabel('Number of demand')
ylabel('Probability')


%% Run Analysis Method
[xx, PI, T] = ExactInDirectPlane(p_sim, p_type, n_sat, Kappa, Q, R, dt_mc, dt_plane);

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
plot(0:length(PI.pi_dmd)-1, PI.pi_dmd, 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

