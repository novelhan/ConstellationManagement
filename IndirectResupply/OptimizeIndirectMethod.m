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

%.. Global Constants
global R_earth mu_earth J2_earth
R_earth = 6400; % [km]
mu_earth = 3.986 * 10^5; % [km^3/sec^2]
J2_earth = 0.00108263;

%.. Conversion Parameters
R2D = 180/pi;   % Radian to Degree
D2R = 1/R2D;

%% System Parameter
%.. Marcov Chain Period
dt_mc = 0.5;                 % [day]

%.. Failure rate
f_ref = 0.1/365; % [#/day]
f_type = 1; % 0:Const, 1:State-Dependant

%.. Analysis Method for Parking Analysis
method = 0; % 0:Time based, 1:Ratio based

%.. Constellation Orbit Parameters (Given)
N_plane     =   40;                     % The number of orbital planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane
h_plane     =   1200;                   % Altitude of in-plane orbits [km]
i_plane     =   50 * D2R;               % Inclination of in-plane orbits [rad]

%.. Parking Orbit Parameters (Given)
i_park      =   50 * D2R;               % Inclination of parking orbits [rad]

%.. Launch Vehicle Lead Time Parameters (Given)
mu_lv_heavy = 60;       % [day]
dt_lv_heavy = 30;       % [day]
N_lv_max_heavy = 40;    % # of max sat. capa. for heavy LV

%.. Launch Vehicle Cost Model
c_build = 0.5; %[$/sat]
c_hold_plane = 5/365; % 0.5 or 5[$/sat/day]
c_hold_park = 5/365; % 0.5 or 5 [$/sat/day]
c_lv_heavy_part = 2; %[$/sat/launch]
c_lv_heavy_full = 67; %[%/launch]
c_transfer = 0.5; % Risk and additional cost for indirect transfer [M$/batch]
c_fuel = 0.001; % Fuel cost for indirect transfer [M$/batch]

%.. Desired System Performance
p_loss_plane  = 0.05; % 0.01 Resilience Prob, the prob. having X < n_sat must be smaller than p_loss
p_loss_park = 0.2; % Resilience Prob, the prob. having X == 0 must be smaller than p_loss

%.. Save the input structure
ParaInPlane.dt_mc = dt_mc;
ParaInPlane.f_ref = f_ref;
ParaInPlane.f_type = f_type;
ParaInPlane.N_sat = N_sat;
ParaInPlane.N_orbit = N_plane;
ParaInPlane.alt = h_plane;
ParaInPlane.inc = i_plane;

ParaParking.dt_mc = dt_mc;
ParaParking.mu_lv = mu_lv_heavy;
ParaParking.dt_lv = dt_lv_heavy;
ParaParking.inc = i_park;
ParaParking.method = method;

ParaCost.c_build = c_build;
ParaCost.c_hold_plane = c_hold_plane;
ParaCost.c_hold_park = c_hold_park;
ParaCost.c_lv_heavy_part = c_lv_heavy_part;
ParaCost.c_lv_heavy_full = c_lv_heavy_full;
ParaCost.c_transfer = c_transfer;
ParaCost.c_fuel = c_fuel;
ParaCost.m_sat = 150; % [kg]
ParaCost.m_bus = 100; % [kg]
ParaCost.Vex = 2.16; % [km/s]

ParaConst.p_loss_plane = p_loss_plane;
ParaConst.p_loss_park = p_loss_park;
ParaConst.N_lv_max_heavy = N_lv_max_heavy;

%% Design Variable
%.. In-Plane Reorder Size Range (Design Var.)
Qi_min = 1;
Qi_max = 10;

%.. In-Plane Reorder Point Range (Design Var.)
Ri_min = N_sat - Qi_max;
Ri_max = N_sat + Qi_max;

%.. Launch Vehicle Size Range (Design Var.)
Qp_min = 1;
Qp_max = N_lv_max_heavy;

%.. Reorder Point Range (Design Var.)
Rp_min = 1;
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
                    'MaxGenerations',1200,'PopulationSize',600,'Display','iter');
rng default % For reproducibility
x_opt = ga(@(x) CostInDirectResupply(x, ParaCost, ParaInPlane, ParaParking), 6, [], [], [], [], LB, UB,...
           @(x) ConstInDirectResupply(x, ParaConst, ParaInPlane, ParaParking), intcon, opts);

%% Check Result
% x_opt = [1.0000   40.0000    9.0000    5.0000    8.0000  574.3653];
% Ref. 4.0000   41.0000   10.0000    2.0000    3.0000  713
% High holding cost. 1.0000   40.0000    9.0000    5.0000    8.0000  574.3653
% Tight const. 2.0000   41.0000   20.0000    2.0000    5.0000  555
% 4.0000   42.0000   10.0000    2.0000    3.0000  571.3087
% 4.0000   41.0000   10.0000    2.0000    4.0000  530.5247
[J, Cost] = CostInDirectResupply(x_opt, ParaCost, ParaInPlane, ParaParking);
c = ConstInDirectResupply(x_opt, ParaConst, ParaInPlane, ParaParking);
disp('External Function Results')
disp(['Total Cost: ',num2str(J)])
disp(['Build Cost: ',num2str(Cost.C_build)])
disp(['Holding Cost: ',num2str(Cost.C_hold)])
disp(['Launch Cost: ',num2str(Cost.C_launch)])
disp(['Transfer Cost: ',num2str(Cost.C_transfer)])
disp(['P(Xi< N_sat): ',num2str(c(1)+p_loss_plane)])
disp(['P(Xp = 0): ',num2str(c(2)+p_loss_park)])
       
