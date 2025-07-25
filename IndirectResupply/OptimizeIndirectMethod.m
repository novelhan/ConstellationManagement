% This script optimizes the indirect spare management strategy using the genetic algorithm.
% The decision variables to optimize are (ri,qi), (rp,qp), Np, hp.
% The cost function is the effective operational cost of the strategy per unit cycle
% The constraints include the payload capacity of LV, resilience of the spare strategy
% (!!!!) The script requires gads_toolbox and optimization_toolbox.

clear all
close all
clc

%.. PATH
mfilepath = pwd;
idcs = strfind(mfilepath,'\');
libdir = mfilepath(1:idcs(end));
addpath([libdir, 'CommonSource'])

%.. Conversion Parameters
R2D = 180/pi;   % Radian to Degree
D2R = 1/R2D;

%% System Parameter
%.. Marcov Chain Period
dt_mc = 0.5;                 % [day]
ParaInPlane.dt_mc = dt_mc;
ParaParking.dt_mc = dt_mc;

%.. Failure rate
% 0.05, 0.1, 0.2
f_fail = 0.1/365; % [#/day]
f_type = 1; % 0:Const, 1:State-Dependant
ParaInPlane.f_sim = f_fail;
ParaInPlane.f_type = f_type;

%.. Constellation Orbit Parameters (Given)
N_plane     =   40;                     % The number of orbital planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane
h_plane     =   1200;                   % Altitude of in-plane orbits [km]
i_plane     =   50 * D2R;               % Inclination of in-plane orbits [rad]

ParaInPlane.n_sat = N_sat;
ParaInPlane.N_orbit = N_plane;
ParaInPlane.alt = h_plane;
ParaInPlane.inc = i_plane;

%.. Parking Orbit Parameters (Given)
i_park      =   50 * D2R;               % Inclination of parking orbits [rad]

%.. Launch Vehicle Lead Time Parameters (Given)
mu_LV   =   60;     % [day]
T0_LV   =   30;     % [day]
ParaParking.mu_lv = mu_LV;
ParaParking.dt_lv = T0_LV;
ParaParking.inc = i_park;
ParaParking.method = 0;

%.. Launch Vehicle Cost Model
c_build = 0.5; %[$/sat]
c_holding = 0.5/365; %[$/sat/day]
c_lv_part = 10; %[$/sat/launch]
lv_full_discount = 0.05; %[%/launch]
Qp_max = 10;
ParaCost.c_build = c_build;
ParaCost.c_holding = c_holding;
ParaCost.c_lv_part = c_lv_part;
ParaCost.lv_full_discount = lv_full_discount;
ParaCost.Qmax = Qp_max;

%.. Desired System Performance
p_loss  =   0.05; % Prob. having less than N_sat number of satellite
ParaConst.p_loss = p_loss;


%% Design Variable
%.. In-Plane Reorder Size Range (Design Var.)
Qi_min = 1;
Qi_max = 10;

%.. In-Plane Reorder Point Range (Design Var.)
Ri_min = N_sat - 5;
Ri_max = N_sat + 5;

%.. Launch Vehicle Size Range (Design Var.)
Qp_min = 2;
Qp_max = Qp_max;

%.. Reorder Point Range (Design Var.)
Rp_min = 0;
Rp_max = Qp_max - 1;

%.. # of Parking Orbit Range
Np_min = 2;
Np_max = 20;

%.. Parking Orbit Altitude Range
hp_min = 500;
hp_max = 1000;


%% Optimize
%.. Lower/Upper Bound
LB = [Qi_min, Ri_min, Qp_min, Rp_min, Np_min, hp_min];
UB = [Qi_max, Ri_max, Qp_max, Rp_max, Np_max, hp_max];

%.. Integer Index
intcon = [1 2 3 4 5];

opts = optimoptions('ga','MaxStallGenerations',50,'FunctionTolerance',1e-10,...
                    'MaxGenerations',500,'PopulationSize',200,'Display','iter');
rng default % For reproducibility
x_opt = ga(@(x) CostFun(x, ParaCost, ParaInPlane, ParaParking), 6, [], [], [], [], LB, UB,...
           @(x) ConstFun(x, ParaConst, ParaInPlane, ParaParking), intcon, opts);

%%
CostFun(x_opt, ParaCost, ParaInPlane, ParaParking)
ConstFun(x_opt, ParaConst, ParaInPlane, ParaParking)

%% Cost Function
function J = CostFun(x, ParaCost, ParaInPlane, ParaParking)
    %.. Set Design Variable
    ParaInPlane.Q = x(1);
    ParaInPlane.R = x(2);
    ParaParking.Q = x(3);
    ParaParking.R = x(4);
    ParaParking.N_orbit = x(5);
    ParaParking.alt = x(6);

    %.. Compute RAAN Drift Period
    [dt_plane, dt_park] = ComputeRaanPeriod(ParaInPlane, ParaParking);
    ParaInPlane.dt_plane = dt_plane;
    ParaParking.dt_park = dt_park;

    %.. Analyze the performance
    [PI_i, PI_p, T_i, T_p] = ExactInDirectProb(100, ParaInPlane, ParaParking);
    
    %.. Cost for making spare satellites per unit time
    Qi = ParaInPlane.Q; 
    Qp = ParaParking.Q; 
    N_plane = ParaInPlane.N_orbit;
    N_park = ParaParking.N_orbit;
    J_build = ParaCost.c_build * N_park * Qp / T_p.T_ir; % c_build: cost per batch
    
    %.. Cost for LV launching
    %.. TODO: [Qmax/Qi] = Qp 인 경우가 full payload가 되도록 변경
    %.. 크기는 신경쓰지 말고, 무게의 총합이 일정 이하가 되도록
    if Qp == ParaCost.Qmax
        J_launch = N_park * (1 - ParaCost.lv_full_discount) * ParaCost.c_lv_part * ParaCost.Qmax / T_p.T_ir;
    else
        J_launch = N_park * ParaCost.c_lv_part * Qp / T_p.T_ir;
    end
    
    %.. Cost for holding redundant spares
    xx_i = length(PI_i.pi_ir)-1:-1:0;
    xx_p = length(PI_p.pi_ir)-1:-1:0;
    idx_i = find(xx_i == ParaInPlane.n_sat + 1);
    idx_p = find(xx_p == 1);
    J_hold_i = N_plane * ParaCost.c_holding * (xx_i(1:idx_i)*PI_i.pi_ir(1:idx_i) - ParaInPlane.n_sat);
    J_hold_p = N_park * ParaCost.c_holding * Qi *(xx_p(1:idx_p)*PI_p.pi_ir(1:idx_p));
    
    %.. Total Cost
    J = J_build + J_launch + J_hold_i + J_hold_p;
end

%% Constraint Function
function [c,ceq] = ConstFun(x, ParaConst, ParaInPlane, ParaParking)
    %.. Set Design Variable
    ParaInPlane.Q = x(1); % Qi
    ParaInPlane.R = x(2); % Ri
    ParaParking.Q = x(3); % Qp
    ParaParking.R = x(4); % Rp
    ParaParking.N_orbit = x(5); % N_park
    ParaParking.alt = x(6); 

    %.. Compute RAAN Drift Period
    [dt_plane, dt_park] = ComputeRaanPeriod(ParaInPlane, ParaParking);
    ParaInPlane.dt_plane = dt_plane;
    ParaParking.dt_park = dt_park;

    %.. Analyze the performance
    [PI_i, PI_p, T_i, T_p] = ExactInDirectProb(100, ParaInPlane, ParaParking);

    %.. p_loss < eps
    xx_i = length(PI_i.pi_ir)-1:-1:0;
    idx = find(xx_i < ParaInPlane.n_sat);
    p_loss = sum(PI_i.pi_ir(idx));
    c = p_loss - ParaConst.p_loss;
    
    %.. Rp + 1 <= Qp
    if ParaParking.method == 0
        c(2) = - ParaParking.Q + ParaParking.R + 1;
    end
    %.. Qi*Qp <= Q_max TODO
    % c(3) = x(1) * x(3) - ParaCost.Q_max;
    
    %.. 0 < Rp < Qp_max (Handled by external feature)
    %.. 

    
    
    ceq = [];
end

%% RAAN Period Computing Function
function [dt_plane, dt_park] = ComputeRaanPeriod(ParaInPlane, ParaParking)
    %.. Earth Parameters
    R_earth     =   6400;
    mu_earth    =   3.986 * 10^5; 
    J2_earth    =   0.00108263;
    
    %.. In-Orbit Parameters
    a_plane     =   R_earth + ParaInPlane.alt;
    e_plane     =   0;
    i_plane     =   ParaInPlane.inc;
    n_plane     =   sqrt(mu_earth/a_plane^3);
    N_plane     =   ParaInPlane.N_orbit;
    Wdot_orbit  =   -3/2 * J2_earth* n_plane * R_earth^2 / (a_plane*(1-e_plane))^2 * cos(i_plane);
    
    %.. Parking-Orbit Parameters
    a_park      =   R_earth + ParaParking.alt;
    e_park      =   0;
    i_park      =   ParaParking.inc;
    n_park      =   sqrt(mu_earth/a_park^3);
    N_park      =   ParaParking.N_orbit;
    Wdot_park   =   -3/2 * J2_earth * n_park * R_earth^2 / (a_park*(1-e_park))^2 * cos(i_park);

    %.. RAAN Drift period for in-plane and parking orbit
    Wdrift      =   abs(Wdot_park - Wdot_orbit);
    dt_plane    =   2*pi/ N_park / Wdrift / 24 / 3600; % RAAN Lead Time [day]
    dt_park     =   2*pi/ N_plane / Wdrift / 24 / 3600; % RAAN Lead Time [day]
end