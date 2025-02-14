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

%% Test Param
% rng('default')
iter_max    =   5;

%.. Sim time
dt_sim      =   1;                  % [day]
time_sim    =   0:dt_sim:365*100;

%.. Conversion Parameters
R2D = 180/pi;   % Radian to Degree
D2R = 1/R2D;

%.. Earth Parameters
R_earth     =   6400;
mu_earth    =   3.986 * 10^5; 
J2_earth    =   0.00108263;

%.. Constellation Parameters
N_plane     =   40;                     % The number of orbital planes for constellation
N_park      =   3;                      % The number of parking planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane
h_plane     =   1200;                   % Altitude of in-plane orbits [km]
h_park      =   600;                    % Altitude of parking orbits [km]
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
Wdrift_dt   =   Wdrift * dt_sim * 24 * 3600; % [rad/dt]

%.. In-Orbit/Parking Period
dt_park     =   2*pi/ N_plane / Wdrift_dt; % Every Contact Period for Parking orbit [dt]
dt_plane    =   2*pi/ N_park / Wdrift_dt; % Every Contact Period for In-Plane orbit [dt]

%.. Indirect LV Parameters (goes to parking orbit)
mu_lv_park      =   60;
dt_lv_park      =   30;                % [day]
cnt_lv_park     =   round(dt_lv_park/dt_sim);

%.. Direct LV Parameters (goes to in-plane orbit)
mu_lv_plane     =   30;
dt_lv_plane     =   15;                % [day]
cnt_lv_plane    =   round(dt_lv_plane/dt_sim);

%.. Failure rate 
%!! If failure rate >> resupply speed -> the method will give poor result
p_fail      =   0.1/365;            % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_type      =   0; % 0 for const, 1 for state dependant

%.. Hybrid In-plane (Q1,R1,Q2,R2) Policy Parameter
Qi1     =   4;
Ri1     =   42;
Qi2     =   2;
Ri2     =   40;

%.. Parking (Q,R) Policy Parameter
kq      =   8;
kr      =   7;
Qp      =   Qi1*kq;
Rp      =   Qi1*kr;

%.. State Parameter
Xi_max = Qi1 + Qi2 + max(Ri1, Ri2); % Max State Level of in-plane
Xi_num = 0:1:Xi_max; % State Counts
Xp_max = kr + kq; % Max State Level of parking
Xp_num = 0:1:Xp_max; % State Counts
Di_max = ceil(Xi_max/Qi1);

% Max Lead Time Bin (mean + 5*sigma) [dt_sim]
LVp_max = ceil((dt_lv_park + mu_lv_park + 5*mu_lv_park)/dt_sim); 
LVi_max = ceil((dt_lv_plane + mu_lv_plane + 5*mu_lv_plane)/dt_sim); 

%% Simulation Result Save
% The number of availalbe stock at each time step / histogram
Ni_on   =   zeros(length(time_sim), N_plane, iter_max);
Np_on   =   zeros(length(time_sim), N_park, iter_max);
Xi_on   =   zeros(Xi_max+1, N_plane, iter_max);
Xp_on   =   zeros(Xp_max+1, N_park, iter_max);

% The histogram of the number of stock right after the resupply moment
Xi_q1   =   zeros(Xi_max+1, N_plane, iter_max);
Xi_q2   =   zeros(Xi_max+1, N_plane, iter_max);
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
%     v_plane2park=   zeros(N_plane, length(time_sim)); % Needed only for simulator vaildation
    for i = 1:N_park
        [min_angle, min_orbit] = min(mod(Wplane0 - Wpark0(i),2*pi));
        T0 = min_angle/Wdrift_dt;   % Time for first RAAN match
        Tf = time_sim(end);
        Tidx = 1 + round((T0 + 0:dt_park:Tf)/dt_sim); % Time index for RAAN match for entire period
        Oidx = mod(min_orbit - 1 + (0:(length(Tidx)-1)), N_plane) + 1; % In-orbit index at Time index
        for j = 1:length(Tidx)
            if time_sim(Tidx(j)) <= Tf
                v_park2plane(i,Tidx(j)) = Oidx(j);
%                 v_plane2park(Oidx(j),Tidx(j)) = i; 
            end  
        end
    end
    
    % Needed only for simulator vaildation 
    if 0
        disp(['Parking Orbit Period should be:' num2str(dt_park*dt_sim)])
        for i = 1:N_park
            disp(['Parking Orbit', num2str(i), ' is ', num2str(mean(diff(find(v_park2plane(1,:)' ~= 0)))*dt_sim)])
        end

        disp(' ')
        disp(['In-Plane Orbit Period should be:' num2str(dt_plane*dt_sim)])
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
                Xi_lv(i,itr) = Xi_lv(i,itr) + 1;

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
                t_Lv = dt_lv_park + CustomExpRnd(mu_lv_park,1);
                LVp_on(j) = ceil( t_Lv/dt_sim );
                
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
                t_Lv = dt_lv_plane + CustomExpRnd(mu_lv_plane,1);
                LVi_on(i) = ceil( t_Lv/dt_sim );
                
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

figure(1)
plot(time_sim, Ni_on(:,:,1))
xlabel('Time [day]')
ylabel('Number of in-plane stock for entire period')

figure(2)
plot(time_sim, Np_on(:,:,1))
xlabel('Time [day]')
ylabel('Number of parking stock for entire period')

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
xlabel('Number of parking stock for entire period')
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
dt_mc       =   dt_sim/2;                 % [day]

%.. InPlane Parameter
ParaFail.dt_mc = dt_mc;
ParaFail.f_ref = p_fail;
ParaFail.f_type = p_type;
ParaFail.n_sat = N_sat;

%.. Direct Resupply Parameter
ParaDirect.mu = mu_lv_plane;
ParaDirect.dt_lv = dt_lv_plane;
ParaDirect.Q2 = Qi2;
ParaDirect.R2 = Ri2;

%.. Indirect Resupply Parameter
ParaIndirect.dt_plane = dt_plane;
ParaIndirect.kappa = ones(Di_max+1,1);
ParaIndirect.Q1 = Qi1;
ParaIndirect.R1 = Ri1;

%.. Full state result
ParaDim.flag = 0;
ParaDim.x_min = min(Ri1,Ri2) - ceil((3*p_fail*N_sat)*(dt_lv_plane + 2*mu_lv_plane));
ParaDim.n_seg = round( (dt_lv_plane/dt_mc - 1)/3 );

%.. Save Input structure for InPlane
ParaInPlane.ParaFail = ParaFail;
ParaInPlane.ParaDirect = ParaDirect;
ParaInPlane.ParaIndirect = ParaIndirect;
ParaInPlane.ParaDim = ParaDim;

%.. Parking Parameter
ParaParking.Q = kq;
ParaParking.R = kr;
ParaParking.dt_mc = dt_mc;
ParaParking.dt_park = dt_park;
ParaParking.dt_lv = dt_lv_park;
ParaParking.mu_lv = mu_lv_park;
ParaParking.method = 1;

[PI_i, PI_p, T_i, T_p, err] = HybridProb(100, ParaInPlane, ParaParking);

pi_hr = mean(PI_i.pi_hr,2);
xx_i = length(pi_hr)-1:-1:0;
xx_p = length(PI_p.pi_ir)-1:-1:0;

figure(3); hold on
plot(xx_i(1:13), pi_hr(1:13), 'r*')
xlim([xx_i(13)-0.5, xx_i(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(4); hold on
plot(xx_p, PI_p.pi_ir, 'r*')
% xlim([xx_p(13)-0.5, xx_p(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(42); hold on
plot(xx_p, PI_p.pi_q, 'r*')
% xlim([xx_p(13)-0.5, xx_p(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(43); hold on
plot(xx_p, PI_p.pi_r, 'r*')
% xlim([xx_p(13)-0.5, xx_p(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(44); hold on
plot(0:length(PI_i.pi_dmd)-1, PI_i.pi_dmd, 'r*')
% xlim([xx_p(13)-0.5, xx_p(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(31)
plot(xx_i(1:13), pi_hr(1:13), 'yo', 'MarkerFaceColor','y')
h = zeros(3, 1);
h(1) = plot(nan,'k*');
h(2) = plot(nan,'ro', 'MarkerFaceColor','r');
h(3) = plot(nan,'yo', 'MarkerFaceColor','y');
legend(h, 'Each Orbit','Avg','Sol');

figure(41)
plot(xx_p, PI_p.pi_ir, 'yo', 'MarkerFaceColor','y')
h = zeros(3, 1);
h(1) = plot(nan,'k*');
h(2) = plot(nan,'ro', 'MarkerFaceColor','r');
h(3) = plot(nan,'yo', 'MarkerFaceColor','y');
legend(h, 'Each Orbit','Avg','Sol');

figure(51); hold on
plot(xx_i(1:13), PI_i.pi_r1(1:13), 'r*')

figure(52); hold on
plot(xx_i(1:13), PI_i.pi_r2(1:13), 'r*')

figure(53); hold on
plot(xx_i(1:13), PI_i.pi_q1(1:13), 'r*')

figure(54); hold on
plot(xx_i(1:13), PI_i.pi_q2(1:13), 'r*')

figure(5)
plot(1:length(err), log10(err), '-o')
xlabel('Iterations')
ylabel('log_{10} of relative error')