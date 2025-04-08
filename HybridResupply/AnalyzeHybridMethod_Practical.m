% This script vaildates the analysis method of hybrid spare management strategy.
% Proposed Markov Chain analysis results are compared with the Monte Carlo simulation results.
% (!!) This script assumes practical configuration for in-plane and parking orbits
% (!!) Test under given Orbital elements (then corresponding Tpark and Tplane are found)

close all
clear all
clc

%.. PATH
mfilepath = pwd;
idcs = strfind(mfilepath,'\');
libdir = mfilepath(1:idcs(end));
addpath([libdir, 'CommonSource'])
addpath([libdir, 'IndirectResupply'])

%.. Global Constants
global R_earth mu_earth J2_earth
R_earth = 6400; % [km]
mu_earth = 3.986 * 10^5; % [km^3/sec^2]
J2_earth = 0.00108263;

% x = [4.0000   40.0000    3.0000   38.0000   10.0000    2.0000    4.0000  600];
% x = [2.0000   41.0000    1.0000   36.0000   20.0000    3.0000    3.0000  530.1923]; % pi_q2 가 이상한 거 확인
x = [1.0000   40.0000    1.0000   37.0000    8.0000    7.0000    7.0000  508.7363];
%% Test Param
% rng('default')
iter_max    =   10;

%.. Sim time
dt_sim      =   0.5;                  % [day]
time_sim    =   0:dt_sim:365*10;

%.. Conversion Parameters
R2D = 180/pi;   % Radian to Degree
D2R = 1/R2D;

%.. Constellation Parameters
N_plane     =   40;                     % The number of orbital planes for constellation
N_park      =   x(7);                      % The number of parking planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane
h_plane     =   1200;                   % Altitude of in-plane orbits [km]
h_park      =   x(8);                    % Altitude of parking orbits [km]
i_plane     =   50 * D2R;               % Inclination of in-plane orbits [rad]
i_park      =   50 * D2R;               % Inclination of parking orbits [rad]

%.. Compute RAAN Drift 
a_plane     =   R_earth + h_plane;
e_plane     =   0;
n_plane     =   sqrt(mu_earth/a_plane^3);
Wdot_plane  =   -3/2 * J2_earth* n_plane * R_earth^2 / (a_plane*(1-e_plane))^2 * cos(i_plane);

a_park      =   R_earth + h_park;
e_park      =   0;
n_park      =   sqrt(mu_earth/a_park^3);
Wdot_park   =   -3/2 * J2_earth * n_park * R_earth^2 / (a_park*(1-e_park))^2 * cos(i_park);

Wdrift      =   abs(Wdot_park - Wdot_plane); % [rad/sec]
Wdrift_dt   =   Wdrift * 24 * 3600; % [rad/day]

%.. In-Orbit/Parking Period
dt_park     =   2*pi/ N_plane / Wdrift_dt; % Every Contact Period for Parking orbit [day]
dt_plane    =   2*pi/ N_park / Wdrift_dt; % Every Contact Period for In-Plane orbit [day]

%.. Indirect LV Parameters (goes to parking orbit)
mu_lv_heavy      =   60;
dt_lv_heavy      =   30;                % [day]
cnt_lv_heavy     =   round(dt_lv_heavy/dt_sim);

%.. Direct LV Parameters (goes to in-plane orbit)
mu_lv_small     =   40;
dt_lv_small     =   20;                % [day]
cnt_lv_small    =   round(dt_lv_small/dt_sim);

%.. Failure rate 
%!! If failure rate >> resupply speed -> the method will give poor result
p_fail      =   0.1/365;            % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_type      =   1; % 0 for const, 1 for state dependant

%.. Hybrid In-plane (Q1,R1,Q2,R2) Policy Parameter
Qi1     =   x(1);
Ri1     =   x(2);
Qi2     =   x(3);
Ri2     =   x(4);

%.. Parking (Q,R) Policy Parameter
kq      =   x(5);
kr      =   x(6);
Qp      =   Qi1*kq;
Rp      =   Qi1*kr;

%.. State Parameter
Xi_max = Qi1 + Qi2 + max(Ri1, Ri2); % Max State Level of in-plane
Xi_num = 0:1:Xi_max; % State Counts
Xp_max = kr + kq; % Max State Level of parking
Xp_num = 0:1:Xp_max; % State Counts
Di_max = ceil(Xi_max/Qi1);

% Max Lead Time Bin (mean + 5*sigma) [dt_sim]
LVp_max = ceil((dt_lv_heavy + mu_lv_heavy + 5*mu_lv_heavy)/dt_sim); 
LVi_max = ceil((dt_lv_small + mu_lv_small + 5*mu_lv_small)/dt_sim); 

%% Simulation Result Save
% The number of availalbe stock at each time step / histogram
Ni_on   =   zeros(length(time_sim), N_plane, iter_max);
Np_on   =   zeros(length(time_sim), N_park, iter_max);
Xi_on   =   zeros(Xi_max+1, N_plane, iter_max);
Xp_on   =   zeros(Xp_max+1, N_park, iter_max);

% The histogram of the number of stock right after the resupply moment
Xi_q1   =   zeros(Xi_max+1, N_plane, iter_max); % Indirect
Xi_q2   =   zeros(Xi_max+1, N_plane, iter_max); % Direct
Xp_q    =   zeros(Xp_max+1, N_park, iter_max);

% The histogram of the number of stock at the reordering moment
Xi_r1   =   zeros(Xi_max+1, N_plane, iter_max);
Xi_r2   =   zeros(Xi_max+1, N_plane, iter_max);
Xp_r    =   zeros(Xp_max+1, N_park, iter_max);

% The histogram of the number of demand stock at the reordering moment
Xi_dmd  =   zeros(Di_max+1, N_plane, iter_max);

% The histogram of the number of N_plane stock
Xi_av = zeros(Di_max+1, Di_max+1, N_plane, iter_max); % dmd, sup, plane, iter
Xp_av = zeros(Di_max+1, Di_max+1, N_park, iter_max); % dmd, sup, park, iter

% The histogram of the number of launch vehicle
Xi_lv = zeros(N_plane, iter_max);
Xp_lv = zeros(N_park, iter_max); 


%% Run Each Simulation
for itr = 1:iter_max
    disp(['Iter: ', num2str(itr)])
    
    % The number of satellite at current time step
    Ni_on_k   =   (Ri1 + round(Qi1*rand(1,N_plane)));
    Np_on_k   =   kr+round(kq*rand(1,N_park));

    % Generate Future Contact Momment
    Wpark0      =   mod(2*pi*rand + 2*pi/N_park*(0:1:N_park-1)',2*pi);
    Wplane0     =   2*pi/N_plane*(0:1:N_plane-1)';
    v_park2plane=   zeros(N_park, length(time_sim));
    v_plane2park=   zeros(N_plane, length(time_sim)); % Needed only for simulator vaildation
    for i = 1:N_park
        [min_angle, min_orbit] = min(mod(Wplane0 - Wpark0(i),2*pi));
        T0 = min_angle/Wdrift_dt;   % Time for first RAAN match
        Tf = time_sim(end);
        Tidx = 1 + round((T0 + 0:dt_park:Tf)/dt_sim); % Time index for RAAN match for entire period
        Oidx = mod(min_orbit - 1 + (0:(length(Tidx)-1)), N_plane) + 1; % In-orbit index at Time index
        for j = 1:length(Tidx)
            if time_sim(Tidx(j)) <= Tf
                v_park2plane(i,Tidx(j)) = Oidx(j);
                v_plane2park(Oidx(j),Tidx(j)) = i; 
            end  
        end
    end
    
    % Needed only for simulator vaildation 
    if 0
        disp(['Parking Orbit Period should be:' num2str(dt_park)])
        for i = 1:N_park
            disp(['Parking Orbit', num2str(i), ' is ', num2str(mean(diff(find(v_park2plane(1,:)' ~= 0)))*dt_sim)])
        end

        disp(' ')
        disp(['In-Plane Orbit Period should be:' num2str(dt_plane)])
        for i = 1:N_plane
            disp(['In-Plane Orbit', num2str(i), ' is ', num2str(mean(diff(find(v_plane2park(1,:)' ~= 0)))*dt_sim)])
        end
    end

    % Set LV Launch Flag (-1 for ready for launch)
    LVi_on = -ones(1,N_plane);
    LVp_on = -ones(1,N_park);

    % Apply Policy
    for k = 1:length(time_sim)
        % Order of simulation
        % In-Plane: Failure -> Check Qi2 arrival -> Apply (Qi1,Ri1) -> Check Ri2
        % Parking: Check Qp arrival -> Demand Distribution -> Check Rp
        
        %%% 1. Generate Sample of In-Plane Failure at Every Time Step
        for i = 1:N_plane
            if p_type == 0 %.. Const Failure Rate
                N_fail = CustomPoisRnd(N_sat*p_sim, 1);
            else %.. State Dependant Failure Rate
                if Ni_on_k(i) > N_sat
                    N_fail = CustomPoisRnd(N_sat*p_sim, 1);
                else
                    N_fail = CustomPoisRnd(Ni_on_k(i)*p_sim, 1);
                end
            end

            % In-Plane Satellites after failure
            Ni_on_k(i) = max([Ni_on_k(i) - N_fail, 0]);
        end

        %%% 2-1. Check Parking Resupply Arrival
        for j = 1:N_park
            if LVp_on(j) == 0 % Resupply has arrived
                % Update Non and Xp_q
                Np_on_k(j) = Np_on_k(j) + kq;
                Xp_q(Np_on_k(j)+1,j,itr) = Xp_q(Np_on_k(j)+1,j,itr) + 1;

                % Update LV Parameters
                LVp_on(j) = -1;

            elseif LVp_on(j) > 0 % Wait for arrival
                LVp_on(j) = LVp_on(j) - 1;
            end
        end

        %%% 2-2. Check In-Plane Direct Resupply Arrival
        for i = 1:N_plane
            if LVi_on(i) == 0 % Resupply has arrived
                % Update Non and Xi_q2
                Ni_on_k(i) = Ni_on_k(i) + Qi2;
                Xi_q2(Ni_on_k(i)+1,i,itr) = Xi_q2(Ni_on_k(i)+1,i,itr) + 1;

                % Update LV Parameters
                LVi_on(i) = -1;

            elseif LVi_on(i) > 0 % Wait for arrival
                LVi_on(i) = LVi_on(i) - 1;
            end
        end
        
        %%% 3. Check In-Orbit Spare Transfer
        for j = 1:N_park
            if v_park2plane(j,k) ~= 0 % RAAN Contact Moment
                % Contacted In-Plane
                i = v_park2plane(j,k);

                % Demand
                n_dmd = max(ceil( (Ri1 + 1 - Ni_on_k(i))/Qi1 ), 0); % Required # of batch for ith in-plane
                n_av = Np_on_k(j); % # of available batch for jth parking
                n_trn = min(n_dmd, n_av); % # of transfered batch for jth park -> ith in-plane

                % Update Xi_r
                Xi_r1(Ni_on_k(i)+1,i,itr) = Xi_r1(Ni_on_k(i)+1,i,itr) + 1;

                % Update Sat
                Ni_on_k(i) = Ni_on_k(i) + Qi1*n_trn;
                Np_on_k(j) = Np_on_k(j) - n_trn;

                % Update Xi_q
                Xi_q1(Ni_on_k(i)+1,i,itr) = Xi_q1(Ni_on_k(i)+1,i,itr) + 1;
                Xi_dmd(n_dmd+1,i,itr) = Xi_dmd(n_dmd+1,i,itr) + 1;
                
                % Update Demand and actual transfer stock
                Xi_av(n_trn+1, n_dmd+1, i, itr) = Xi_av(n_trn+1, n_dmd+1, i, itr) + 1; 
                Xp_av(n_trn+1, n_dmd+1, j, itr) = Xp_av(n_trn+1, n_dmd+1, j, itr) + 1;
            end
        end

        %%% 4-1. Check Parking Reorder
        for j = 1:N_park
            % LV must be available, # of spares < kr 
            if LVp_on(j) == -1 && Np_on_k(j) <= kr
                % Sample Lead Time and Save
                t_LV = dt_lv_heavy + CustomExpRnd(mu_lv_heavy,1);
                LVp_on(j) = ceil( t_LV/dt_sim );
                
                % Update Xp_lv
                Xp_lv(j,itr) = Xp_lv(j,itr) + 1;

                % Update Xr
                Xp_r(Np_on_k(j)+1,j) = Xp_r(Np_on_k(j)+1,j) + 1; % +1 for index
                
                % Save Remaining time step before arrival
                LVp_on(j) = LVp_on(j) - 1;
            end
        end

        %%% 4-2. Check In-plane Direct Reorder
        for i = 1:N_plane
            % LVi must be available, # of spares < Ri2 
            if LVi_on(i) == -1 && Ni_on_k(i) <= Ri2
                % Sample Lead Time and Save
                t_LV = dt_lv_small + CustomExpRnd(mu_lv_small,1);
                LVi_on(i) = ceil( t_LV/dt_sim );
                
                % Update Xp_lv
                Xi_lv(i,itr) = Xi_lv(i,itr) + 1;

                % Update Xi_r2
                Xi_r2(Ni_on_k(i)+1,i) = Xi_r2(Ni_on_k(i)+1,i) + 1; % +1 for index
                
                % Save Remaining time step before arrival
                LVi_on(i) = LVi_on(i) - 1;
            end
        end
        
        %%% 5. Save Stock Profile
        Ni_on(k,:,itr) = Ni_on_k;
        Np_on(k,:,itr) = Np_on_k;
        for i = 1:N_plane
            Xi_on(Ni_on_k(i)+1,i,itr) = Xi_on(Ni_on_k(i)+1,i,itr) + 1;
        end
        for j = 1:N_park
            Xp_on(Np_on_k(j)+1,j,itr) = Xp_on(Np_on_k(j)+1,j,itr) + 1;
        end
        
    end
end


%% Plot Simulation Result
close all

xxi_edge = -0.5:1:(Xi_max+0.5);
xxp_edge = -0.5:1:(Xp_max+0.5);
dmd_edge = -0.5:1:(Di_max+0.5);

%.. Average Result for each different iteration
Xi_on = sum(Xi_on, 3);
Xp_on = sum(Xp_on, 3);
Xp_q = sum(Xp_q,3);
Xp_r = sum(Xp_r,3);
Xi_dmd = sum(Xi_dmd,3);
Xi_r1 = sum(Xi_r1,3);
Xi_r2 = sum(Xi_r2,3);
Xi_q1 = sum(Xi_q1,3);
Xi_q2 = sum(Xi_q2,3);
Xi_av = sum(Xi_av,4);
Xp_av = sum(Xp_av,4);

% figure(1)
% plot(time_sim, Ni_on(:,:,1))
% xlabel('Time [day]')
% ylabel('Number of in-plane stock for entire period')
% 
% figure(2)
% plot(time_sim, Np_on(:,:,1))
% xlabel('Time [day]')
% ylabel('Number of parking stock for entire period')

figure(3); hold on
histogram('BinEdges', xxi_edge, 'BinCounts', sum(Xi_on,2),'Normalization','probability')
xlabel('Number of in-plane stock for entire period')
ylabel('Probability')

Pi_on = zeros(Xi_max+1, N_plane);
for i = 1:N_plane
    Pi_on(:,i) = Xi_on(:,i)/sum(Xi_on(:,i));
end

figure(31); hold on
for i = Ri1-5:Xi_max
    plot(i,Pi_on(i+1,:),'k*')
    plot(i,mean(Pi_on(i+1,:)),'ro','MarkerFaceColor','r')
end
xlabel('Number of in-plane stock for entire period')
ylabel('Probability')

figure(4); hold on
histogram('BinEdges', xxp_edge, 'BinCounts', sum(Xp_on,2),'Normalization','probability')
xlabel('Number of parking stock for entire period')
ylabel('Probability')

Pp_on = zeros(Xp_max+1, N_park);
for i = 1:N_park
    Pp_on(:,i) = Xp_on(:,i)/sum(Xp_on(:,i));
end

figure(41); hold on
for i = 0:Xp_max
    plot(i,Pp_on(i+1,:),'k*')
    plot(i,mean(Pp_on(i+1,:)),'ro')
end
xlabel('Number of parking stock for entire period')
ylabel('Probability')

figure(42); hold on
histogram('BinEdges', xxp_edge, 'BinCounts', sum(Xp_q,2),'Normalization','probability')
xlabel('Number of parking stock right after replenishment')
ylabel('Probability')

figure(43); hold on
histogram('BinEdges', xxp_edge, 'BinCounts', sum(Xp_r,2),'Normalization','probability')
xlabel('Number of parking stock at the reorder moment')
ylabel('Probability')

figure(44); hold on
histogram('BinEdges', dmd_edge, 'BinCounts', sum(Xi_dmd,2),'Normalization','probability')
xlabel('Demand Probability Distribution')
ylabel('Probability')

figure(51); hold on
histogram('BinEdges', xxi_edge, 'BinCounts', sum(Xi_r1,2),'Normalization','probability')
xlabel('Number of in-plane stock at the R1 reorder moment')
ylabel('Probability')

figure(52); hold on
histogram('BinEdges', xxi_edge, 'BinCounts', sum(Xi_r2,2),'Normalization','probability')
xlabel('Number of in-plane stock at the R2 reorder moment')
ylabel('Probability')

figure(53); hold on
histogram('BinEdges', xxi_edge, 'BinCounts', sum(Xi_q1,2),'Normalization','probability')
xlabel('Number of in-plane stock right after Q1 replenishment')
ylabel('Probability')

figure(54); hold on
histogram('BinEdges', xxi_edge, 'BinCounts', sum(Xi_q2,2),'Normalization','probability')
xlabel('Number of in-plane stock right after Q2 replenishment')
ylabel('Probability')

%% Run Analysis Method
%.. Marcov Chain Period
dt_mc       =   dt_sim;                 % [day]

%.. InPlane Parameter
ParaInPlane.dt_mc = dt_mc;
ParaInPlane.f_ref = p_fail;
ParaInPlane.f_type = p_type;
ParaInPlane.N_sat = N_sat;

%.. Direct Resupply Parameter
ParaInPlane.mu = mu_lv_small;
ParaInPlane.dt_lv = dt_lv_small;
ParaInPlane.Q2 = Qi2;
ParaInPlane.R2 = Ri2;

%.. Indirect Resupply Parameter
ParaInPlane.dt_plane = dt_plane;
ParaInPlane.Q1 = Qi1;
ParaInPlane.R1 = Ri1;

%.. Full state result
ParaInPlane.dim_flag = 0;
% ParaInPlane.x_min = min(Ri1,Ri2) - ceil((3*p_fail*N_sat)*(dt_lv_small + 2*mu_lv_small));
% ParaInPlane.n_seg = round( (dt_lv_small/dt_mc - 1)/3 );

%.. Parking Parameter
ParaParking.Q = kq;
ParaParking.R = kr;
ParaParking.dt_mc = dt_mc;
ParaParking.dt_park = dt_park;
ParaParking.dt_lv = dt_lv_heavy;
ParaParking.mu_lv = mu_lv_heavy;
ParaParking.method = 1;

[PI_i, PI_p, T_i, T_p, err] = SolveHybridProb(100, ParaInPlane, ParaParking);

figure(3); hold on
plot(PI_i.X, abs(PI_i.pi_avg), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(4); hold on
plot(PI_p.X, abs(PI_p.pi_avg), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(42); hold on
plot(PI_p.X, abs(PI_p.pi_q), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(43); hold on
plot(PI_p.X, abs(PI_p.pi_r), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(44); hold on
plot(0:length(PI_i.pi_dmd)-1, abs(PI_i.pi_dmd), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(31)
plot(PI_i.X, PI_i.pi_avg, 'yo', 'MarkerFaceColor','y')
h = zeros(3, 1);
h(1) = plot(nan,'k*');
h(2) = plot(nan,'ro', 'MarkerFaceColor','r');
h(3) = plot(nan,'yo', 'MarkerFaceColor','y');
legend(h, 'Each Orbit','Avg','Sol');

figure(41)
plot(PI_p.X, PI_p.pi_avg, 'yo', 'MarkerFaceColor','y')
h = zeros(3, 1);
h(1) = plot(nan,'k*');
h(2) = plot(nan,'ro', 'MarkerFaceColor','r');
h(3) = plot(nan,'yo', 'MarkerFaceColor','y');
legend(h, 'Each Orbit','Avg','Sol');

figure(51); hold on
plot(PI_i.X, abs(PI_i.pi_r1), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(52); hold on
plot(PI_i.X, abs(PI_i.pi_r2), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(53); hold on
plot(PI_i.X, abs(PI_i.pi_q1), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

figure(54); hold on
plot(PI_i.X, abs(PI_i.pi_q2), 'r*')
legend('Sim.', 'Sol.', 'location', 'best')

% figure(5)
% plot(1:length(err), log10(err), '-o')
% xlabel('Iterations')
% ylabel('log_{10} of relative error')


%% Additional Cost Validation
c_build = 0.5; %[$/sat]
c_hold_plane = .5/365; %[$/sat/day]
c_hold_park = .5/365; %[$/sat/day]
c_lv_heavy_part = 2; %[$/sat/launch]
c_lv_heavy_full = 67; %[%/launch]
c_lv_small_part = 3; %[$/sat/launch]
c_lv_small_full = 7.5; %[%/launch]
c_transfer = 0.5; % Risk and additional cost for indirect transfer [M$/batch]
c_fuel = 0.001; % Fuel cost for indirect transfer [M$/batch]
N_lv_max_heavy = 40;
N_lv_max_small = 3;

ParaCost.c_build = c_build;
ParaCost.c_hold_plane = c_hold_plane;
ParaCost.c_hold_park = c_hold_park;
ParaCost.c_lv_heavy_part = c_lv_heavy_part;
ParaCost.c_lv_heavy_full = c_lv_heavy_full;
ParaCost.c_lv_small_part = c_lv_small_part;
ParaCost.c_lv_small_full = c_lv_small_full;
ParaCost.c_transfer = c_transfer;
ParaCost.c_fuel = c_fuel;
ParaCost.m_sat = 150; % [kg]
ParaCost.m_bus = 100; % [kg]
ParaCost.Vex = 2.16; % [km/s]

ParaConst.p_loss_plane = 0;
ParaConst.p_loss_park = 0;
ParaConst.N_lv_max_heavy = N_lv_max_heavy;
ParaConst.N_lv_max_small = N_lv_max_small;

ParaInPlane.alt = h_plane;
ParaInPlane.inc = i_plane;
ParaInPlane.N_orbit = N_plane;
ParaInPlane.dim_flag = 1;

ParaParking.mu_lv = mu_lv_heavy;
ParaParking.dt_lv = dt_lv_heavy;
ParaParking.inc = i_park;
ParaParking.method = 1;

%.. Simulation Result
Si_on = Ni_on - N_sat;
Si_on(Si_on < 0) = 0; 
Ni_avg = mean(mean(mean(Ni_on,3),1),2);
Si_avg = mean(mean(mean(Si_on,3),1),2);
Sp_avg = mean(mean(mean(Np_on,3),1),2);
Np_lv = sum(sum(Xp_lv,2));
Ni_lv = sum(sum(Xi_lv,2));
Ntrn = sum(sum(Xi_av,3),2);

fp_lv = Np_lv/(iter_max*time_sim(end));
fi_lv = Ni_lv/(iter_max*time_sim(end));
f_sat = N_park*Qp*fp_lv + N_plane*Qi2*fi_lv;
f_trn = dot(0:Di_max,Ntrn)/(iter_max*time_sim(end));
p_loss = sum(Xi_on(1:N_sat,:))./sum(Xi_on,1);

disp('Simulation Results')
disp(['Avg. # InPlane Sat: ',num2str(Ni_avg)])
disp(['Avg. # InPlane Spares: ',num2str(Si_avg)])
disp(['Avg. # Parking Spares: ',num2str(Sp_avg)])
disp(['# of bulit sat per unit time: ', num2str(f_sat)])
disp(['# of Small LV per unit time: ', num2str(fi_lv)])
disp(['# of Heavy LV per unit time: ', num2str(fp_lv)])
disp(['# of transfered spares per unit time: ', num2str(f_trn)])
disp(['P(Xi< N_sat): ', num2str(mean(p_loss))])

%.. Analysis result
ni_avg = sum(PI_i.X.*PI_i.pi_avg);
idx = find(PI_i.X == N_sat);
ni_spare = sum((PI_i.X(1:idx-1)-N_sat).*PI_i.pi_avg(1:idx-1));
np_avg = sum(PI_p.X.*PI_p.pi_avg);
fp_lv = N_park/T_p.T_avg;
fi_lv = N_plane/T_i.T_q2;
f_sat = N_park*Qp*fp_lv + N_plane*Qi2*fi_lv;
f_trn = N_plane*(dot(PI_i.X, PI_i.pi_q1) - dot(PI_i.X, PI_i.pi_r1))/Qi1/T_i.T_q1;
p_loss = sum(PI_i.pi_avg(idx+1:end));

disp(' ')
disp('Analysis Results')
disp(['Avg. # Sat: ', num2str(ni_avg)])
disp(['Avg. # InPlane Spares: ',num2str(ni_spare)])
disp(['Avg. # Parking Spares: ',num2str(np_avg)])
disp(['# of bulit sat per unit time: ', num2str(f_sat)])
disp(['# of Small LV per unit time: ', num2str(fi_lv)])
disp(['# of Heavy LV per unit time: ', num2str(fp_lv)])
disp(['# of transfered spares per unit time: ', num2str(f_trn)])
disp(['P(Xi< N_sat): ', num2str(p_loss)])

%.. Using external function
[J, Cost] = CostHybridResupply([Qi1 Ri1 Qi2 Ri2 kq kr N_park h_park], ParaCost, ParaInPlane, ParaParking);
c = ConstHybridResupply([Qi1 Ri1 Qi2 Ri2 kq kr N_park h_park], ParaConst, ParaInPlane, ParaParking);

disp(' ')
disp('External Function Results')
disp(['Total Cost: ',num2str(J)])
disp(['Build Cost: ',num2str(Cost.C_build)])
disp(['Holding Cost: ',num2str(Cost.C_hold)])
disp(['Launch Cost: ',num2str(Cost.C_launch)])
disp(['Transfer Cost: ',num2str(Cost.C_transfer)])
disp(['P(Xi< N_sat): ',num2str(c(1))])
disp(['P(Xp = 0): ',num2str(c(2))])