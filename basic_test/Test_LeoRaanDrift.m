close all
clear all
clc

% This script tests the orbital tranfer feasibility check prediction between parking orbit and
% constellation orbit.

%.. Conversion Parameters
Y2D = 365;      % Year to Day
D2Y = 1/Y2D;
R2D = 180/pi;   % Radian to Degree
D2R = 1/R2D;

%.. Sim Time
T_f         =   10*Y2D; % Simulation Period [Day]
dt_sim      =   0.5; % Time step for simulation [Day]
T_sim       =   0:dt_sim:T_f;

%.. Earth Parameters
R_earth     =   6400;
mu_earth    =   3.986 * 10^5; 
J2_earth    =   0.00108263;

%.. In-Orbit Parameters
a_orbit     =   R_earth + 1200;
e_orbit     =   0;
i_orbit     =   50 * R2D;
n_orbit     =   sqrt(mu_earth/a_orbit^3);
N_orbit     =   40;
Wdot_orbit  =   -3/2 * J2_earth* n_orbit * R_earth^2 / (a_orbit*(1-e_orbit))^2 * cos(i_orbit);

%.. Parking-Orbit Parameters
a_park      =   R_earth + 500;
e_park      =   0;
i_park      =   50 * R2D;
n_park      =   sqrt(mu_earth/a_park^3);
N_park      =   3;
Wdot_park   =   -3/2 * J2_earth * n_park * R_earth^2 / (a_park*(1-e_park))^2 * cos(i_park);

%.. Initial Parking/In Orbit Distribution
Wpark0      =   2*pi/N_park*(0:1:N_park-1)' + 2*pi*rand;
Wpark0      =   wrapTo2Pi(Wpark0);
Worbit0     =   2*pi/N_orbit*(0:1:N_orbit-1)';

%.. Relative RAAN Drift
Wdrift      =   abs(Wdot_park - Wdot_orbit); % [rad/sec]
Wdrift_dt   =   Wdrift * dt_sim * 24 * 3600; % [rad/dt]

%.. Time for parking orbit to rotate i th in orbit to i+1 th in orbit
dt_park     =   2*pi/ N_orbit / Wdrift_dt; % Every Contact Period for Parking orbit [dt]
dt_orbit    =   2*pi/ N_park / Wdrift_dt; % Every Contact Period for In-Plane orbit [dt]

v_i2j = zeros(N_park, length(T_sim)); % Index matrix of i th parking orbit at j th time index (Value = in orbit index)
v_j2i = zeros(N_orbit, length(T_sim));
for i = 1:N_park
    [min_angle, min_orbit] = min(mod(Worbit0 - Wpark0(i),2*pi));
    T0 = min_angle/Wdrift_dt;   % Time for first RAAN match
    Tidx = 1 + round(T0 + 0:dt_park:(T_f/dt_sim)); % Time index for RAAN match for entire period
    Oidx = mod(min_orbit - 1 + (0:(length(Tidx)-1)), N_orbit) + 1; % In-orbit index at Time index
    for j = 1:length(Tidx)
        if Tidx(j) <= (T_f/dt_sim)
            v_i2j(i,Tidx(j)) = Oidx(j);
            v_j2i(Oidx(j),Tidx(j)) = i; 
        end        
    end
end

disp(['Parking Orbit Period should be:' num2str(dt_park*dt_sim)])
for i = 1:N_park
    disp(['Parking Orbit', num2str(i), ' is ', num2str(mean(diff(find(v_i2j(1,:)' ~= 0)))*dt_sim)])
end

disp(' ')
disp(['In-Plane Orbit Period should be:' num2str(dt_orbit*dt_sim)])
for i = 1:N_orbit
    disp(['In-Plane Orbit', num2str(i), ' is ', num2str(mean(diff(find(v_j2i(1,:)' ~= 0)))*dt_sim)])
end

figure(1)
plot(T_sim, v_i2j(1,:),'o')

figure(2)
plot(T_sim, v_i2j(round(N_park/2),:),'o')

figure(3)
plot(T_sim, v_i2j(end,:),'o')
