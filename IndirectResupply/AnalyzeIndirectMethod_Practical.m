% This script vaildates the analysis method of indirect spare management strategy.
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

x = [10.0000   39.0000    2.0000         0    6.0000  500.3649];
%% Test Param
% rng('default')
iter_max    =   2;

%.. Sim time
dt_sim      =   1;                  % [day]
time_sim    =   0:dt_sim:365*100;

%.. Marcov Chain Period
dt_mc       =   dt_sim/2;                 % [day]

%.. Conversion Parameters
R2D = 180/pi;   % Radian to Degree
D2R = 1/R2D;

%.. Earth Parameters
R_earth     =   6400;
mu_earth    =   3.986 * 10^5; 
J2_earth    =   0.00108263;

%.. Constellation Parameters
N_plane     =   40;                     % The number of orbital planes for constellation
N_park      =   x(5);                      % The number of parking planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane
h_plane     =   1200;                   % Altitude of in-plane orbits [km]
h_park      =   x(6);                    % Altitude of parking orbits [km]
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

%.. LV Lead Time
mu_lv       =   60;
dt_lv       =   30;                % [day]
cnt_lv      =   round(dt_lv/dt_sim);

%.. Failure rate 
%!! If failure rate >> resupply speed -> the method will give poor result
p_fail      =   0.1/365;            % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_type      =   0; % 0 for const, 1 for state dependant

%.. In-Plane (Q,R) Policy Parameter
Qi      =   x(1);
Ri      =   x(2);

%.. Parking (Q,R) Policy Parameter
kq      =   x(3);
kr      =   x(4);
Qp      =   Qi*kq;
Rp      =   Qi*kr;

%.. State Parameter
Xi_max = Qi + Ri; % Max State Level of in-plane
Xi_num = 0:1:Xi_max; % State Counts
Xp_max = kr + kq; % Max State Level of parking
Xp_num = 0:1:Xp_max; % State Counts
Di_max = ceil(Xi_max/Qi);
Lv_max = ceil((dt_lv + mu_lv + 5*mu_lv)/dt_sim); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]

%% Simulation Result Save
% The number of availalbe stock at each time step / histogram
Ni_on   =   zeros(length(time_sim), N_plane, iter_max);
Np_on   =   zeros(length(time_sim), N_park, iter_max);
Xi_on   =   zeros(Xi_max+1, N_plane, iter_max);
Xp_on   =   zeros(Xp_max+1, N_park, iter_max);

% The histogram of the number of stock right after the resupply moment
Xi_q    =   zeros(Xi_max+1, N_plane, iter_max);
Xp_q    =   zeros(Xp_max+1, N_park, iter_max);

% The histogram of the number of stock at the reordering moment
Xi_r    =   zeros(Xi_max+1, N_plane, iter_max);
Xp_r    =   zeros(Xp_max+1, N_park, iter_max);

% The histogram of the number of demand stock at the reordering moment
Xi_dmd  =   zeros(Di_max+1, N_plane, iter_max);

% The histogram of the number of available stock
Xi_av = zeros(Di_max+1, Di_max+1, N_plane, iter_max); % dmd, sup, plane, iter
Xp_av = zeros(Di_max+1, Di_max+1, N_park, iter_max); % dmd, sup, park, iter

% The histogram of the number of launch vehicle
Xp_lv = zeros(N_park, iter_max); 


%% Run Each Simulation
for itr = 1:iter_max
    % The number of satellite at current time step
    Ni_on_k   =   (Ri + round(Qi*rand(1,N_plane)));
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
    Lv_on = -ones(1,N_park);

    % Apply Policy
    for k = 1:length(time_sim)

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

        %%% 2. Check Parking Resupply Arrival
        for j = 1:N_park
            if Lv_on(j) == 0 % Resupply has arrived
                % Update Non and Xq
                Np_on_k(j) = Np_on_k(j) + kq;
                Xp_q(Np_on_k(j)+1,j,itr) = Xp_q(Np_on_k(j)+1,j,itr) + 1;

                % Update LV Parameters
                Lv_on(j) = -1;
                Xp_lv(j,itr) = Xp_lv(j,itr) + 1;

            elseif Lv_on(j) > 0 % Wait for arrival
                Lv_on(j) = Lv_on(j) - 1;
            end
        end

        %%% 3. Check In-Orbit Spare Transfer
        for j = 1:N_park
            if v_park2plane(j,k) ~= 0 % RAAN Contact Moment
                % Contacted In-Plane
                i = v_park2plane(j,k);

                % Demand
                n_dmd = max(ceil( (Ri + 1 - Ni_on_k(i))/Qi ), 0); % Required # of batch for ith in-plane
                n_av = Np_on_k(j); % # of available batch for jth parking
                n_trn = min(n_dmd, n_av); % # of transfered batch for jth park -> ith in-plane

                % Update Xi_r
                Xi_r(Ni_on_k(i)+1,i,itr) = Xi_r(Ni_on_k(i)+1,i,itr) + 1;

                % Update Sat
                Ni_on_k(i) = Ni_on_k(i) + Qi*n_trn;
                Np_on_k(j) = Np_on_k(j) - n_trn;

                % Update Xi_q
                Xi_q(Ni_on_k(i)+1,i,itr) = Xi_q(Ni_on_k(i)+1,i,itr) + 1;
                Xi_dmd(n_dmd+1,i,itr) = Xi_dmd(n_dmd+1,i,itr) + 1;
                
                % Update Demand and actual transfer stock
                Xi_av(n_trn+1, n_dmd+1, i, itr) = Xi_av(n_trn+1, n_dmd+1, i, itr) + 1; 
                Xp_av(n_trn+1, n_dmd+1, j, itr) = Xp_av(n_trn+1, n_dmd+1, j, itr) + 1;
            end
        end

        %%% 4. Check Parking Reorder
        for j = 1:N_park
%             LV must be available, RAAN Contact Moment, # of spares < kr 
            if Lv_on(j) == -1 && v_park2plane(j,k) ~= 0 && Np_on_k(j) <= kr % Will save time (If we use Q >= R)
%             if Lv_on(j) == -1 && Np_on_k(j) <= kr
                % Sample Lead Time and Save
                t_Lv = dt_lv + CustomExpRnd(mu_lv,1);
                Lv_on(j) = ceil( t_Lv/dt_sim );
                
                % Update Xp_lv
                Xp_lv(j,itr) = Xp_lv(j,itr) + 1;

                % Update Xr
                Xp_r(Np_on_k(j)+1,j,itr) = Xp_r(Np_on_k(j)+1,j,itr) + 1; % +1 for index
                
                % Save Remaining time step before arrival
                Lv_on(j) = Lv_on(j) - 1;
            end
        end

        % Save Stock Profile
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
Xi_on = sum(Xi_on, 3);
Xp_on = sum(Xp_on, 3);
Xp_q = sum(Xp_q,3);
Xp_r = sum(Xp_r,3);
Xi_dmd = sum(Xi_dmd,3);

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
for i = Ri-5:Xi_max
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
xlabel('Number of batch demanded from in-plane orbits')
ylabel('Probability')


%% Run Analysis Method
%.. InPlane Parameter
ParaInPlane.Q = Qi;
ParaInPlane.R = Ri;
ParaInPlane.n_sat = N_sat;
ParaInPlane.dt_mc = dt_mc;
ParaInPlane.dt_plane = dt_plane;
ParaInPlane.f_sim = p_fail;
ParaInPlane.f_type = p_type;

%.. Parking Parameter
ParaParking.Q = kq;
ParaParking.R = kr;
ParaParking.dt_mc = dt_mc;
ParaParking.dt_park = dt_park;
ParaParking.dt_lv = dt_lv;
ParaParking.mu_lv = mu_lv;
ParaParking.method = 0;

[PI_i, PI_p, T_i, T_p, err] = ExactInDirectProb(100, ParaInPlane, ParaParking);

xx_i = length(PI_i.pi_ir)-1:-1:0;
xx_p = length(PI_p.pi_ir)-1:-1:0;

figure(3); hold on
plot(xx_i(1:13), PI_i.pi_ir(1:13), 'r*')
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
plot(0:Di_max-1, PI_i.pi_dmd, 'r*')
% xlim([xx_p(13)-0.5, xx_p(1)+0.5])
legend('Sim.', 'Sol.', 'location', 'best')

figure(31)
plot(xx_i(1:10), PI_i.pi_ir(1:10), 'yo', 'MarkerFaceColor','y')
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

figure(5)
plot(1:length(err), log10(err), '-o')
xlabel('Iterations')
ylabel('log_{10} of relative error')