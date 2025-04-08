% This script optimizes the hybrid spare management strategy using the genetic algorithm.
% The decision variables to optimize are (ri1,qi1), (ri2,qi2), (rp,qp), Np, hp.
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
addpath([libdir, 'IndirectResupply'])

%.. Conversion Parameters
R2D = 180/pi;   % Radian to Degree
D2R = 1/R2D;

%.. Global Constants
global R_earth mu_earth J2_earth
R_earth = 6400; % [km]
mu_earth = 3.986 * 10^5; % [km^3/sec^2]
J2_earth = 0.00108263;

%% System Parameter
%.. Marcov Chain Period
dt_mc = 0.5;                 % [day]

%.. Failure rate
f_ref = 0.1/365; % [#/day]
f_type = 1; % 0:Constant failure, 1:State Dependent failure

%.. Analysis Method for Parking Analysis
method = 0; % 0:Time based, 1:Ratio based

%.. Constellation Orbit Parameters (Given)
N_plane     =   40;                     % The number of orbital planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane
h_plane     =   1200;                   % Altitude of in-plane orbits [km]
i_plane     =   50 * D2R;               % Inclination of in-plane orbits [rad]

%.. Parking Orbit Parameters (Given)
i_park      =   50 * D2R;               % Inclination of parking orbits [rad]

%.. Indirect LV Parameters (goes to parking orbit)
mu_lv_heavy      =   60;
dt_lv_heavy      =   30;                % [day]
N_lv_max_heavy = 40; %[#sat] max capacitiy for indirect LV

%.. Direct LV Parameters (goes to in-plane orbit)
mu_lv_small     =   40;
dt_lv_small     =   20;                % [day]
N_lv_max_small = 3; %[#sat] max capacitiy for direct LV

%.. Cost Model
c_build = 0.5; % Manufacturing cost per satellite [M$/sat]
c_hold_plane = 5/365; % Holding cost per satellite in In-Plane orbit [M$/sat/day]
c_hold_park = 5/365; %Holding cost per satellite in Parking orbit [M$/sat/day] [M$/sat/day]
c_transfer = 0.5; % Risk and additional cost for indirect transfer [M$/batch]
c_fuel = 0.001; % Fuel cost for indirect transfer [M$/batch]

c_lv_small_part = 3; %[M$/sat/launch] for direct LV
c_lv_small_full = 7.5;%[M$/launch] for direct LV
c_lv_heavy_part = 2; %[M$/sat/launch] for indirect LV
c_lv_heavy_full = 67; %[M$/launch] for indirect LV

%.. Desired System Performance
p_loss_plane  = 0.05; % 0.01 Resilience Prob, the prob. having X < n_sat must be smaller than p_loss
p_loss_park = 0.2; % Resilience Prob, the prob. having X == 0 must be smaller than p_loss

%.. Save the input structure
ParaInPlane.dt_mc = dt_mc;
ParaInPlane.f_ref = f_ref;
ParaInPlane.f_type = f_type;
ParaInPlane.N_sat = N_sat;
ParaInPlane.N_orbit = N_plane;
ParaInPlane.mu = mu_lv_small;
ParaInPlane.dt_lv = dt_lv_small;
ParaInPlane.alt = h_plane;
ParaInPlane.inc = i_plane;
ParaInPlane.dim_flag = 1;
% ParaInPlane.x_min = min(Ri1,Ri2) - ceil((3*p_fail*N_sat)*(dt_lv_small + 2*mu_lv_small));
% ParaInPlane.n_seg = round( (dt_lv_small/dt_mc - 1)/3 );

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
ParaCost.c_lv_small_part = c_lv_small_part;
ParaCost.c_lv_small_full = c_lv_small_full;
ParaCost.c_transfer = c_transfer;
ParaCost.c_fuel = c_fuel;
ParaCost.m_sat = 150; % [kg]
ParaCost.m_bus = 100; % [kg]
ParaCost.Vex = 2.16; % [km/s]

ParaConst.p_loss_plane = p_loss_plane;
ParaConst.p_loss_park = p_loss_park;
ParaConst.N_lv_max_heavy = N_lv_max_heavy;
ParaConst.N_lv_max_small = N_lv_max_small;

%% Design Variable
%.. In-Plane InDirect Reorder Size Range (Design Var.)
Qi1_min = 1;
Qi1_max = 10;

%.. In-Plane InDirect Reorder Point Range (Design Var.)
Ri1_min = N_sat - 10;
Ri1_max = N_sat + 5;

%.. In-Plane Direct Reorder Size Range (Design Var.)
Qi2_min = 1;
Qi2_max = N_lv_max_small;

%.. In-Plane Direct Reorder Point Range (Design Var.)
Ri2_min = N_sat - 10;
Ri2_max = N_sat + 5;

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
LB = [Qi1_min, Ri1_min, Qi2_min, Ri2_min, Qp_min, Rp_min, Np_min, hp_min];
UB = [Qi1_max, Ri1_max, Qi2_max, Ri2_max, Qp_max, Rp_max, Np_max, hp_max];

%.. Integer Index and Opt Option
intcon = [1 2 3 4 5 6 7];
opts = optimoptions('ga','MaxStallGenerations',50,'FunctionTolerance',1e-10,...
                    'MaxGenerations',2000,'PopulationSize',800,'Display','iter');

%.. Run Opt
% rng default % For reproducibility
% x_opt = ga(@(x) CostHybridResupply(x, ParaCost, ParaInPlane, ParaParking), 8, [], [], [], [], LB, UB,...
%            @(x) ConstHybridResupply(x, ParaConst, ParaInPlane, ParaParking), intcon, opts);

x_opt = [ 1.0000   40.0000    1.0000   37.0000    8.0000    7.0000    7.0000  508.7363]; 
% 2.0000   41.0000    1.0000   36.0000   20.0000    3.0000    3.0000  530.1923 Ref.
% 4.0000   41.0000    2.0000   31.0000   10.0000    2.0000    4.0000  530.5247 e_tol
% 4.0000   42.0000    2.0000   32.0000   10.0000    2.0000    3.0000  571.1405
% 1.0000   40.0000    1.0000   30.0000    7.0000    5.0000    9.0000  639.4550 High Holding
% 1.0000   39.0000    2.0000   40.0000    2.0000    1.0000    8.0000  594.3996
% 1.0000   40.0000    3.0000   40.0000    5.0000    2.0000    2.0000  713.9099
% 1.0000   40.0000    1.0000   37.0000    8.0000    7.0000    7.0000  508.7363
%% Check Result
[J, Cost] = CostHybridResupply(x_opt, ParaCost, ParaInPlane, ParaParking);
c = ConstHybridResupply(x_opt, ParaConst, ParaInPlane, ParaParking);

disp('External Function Results')
disp(['Total Cost: ',num2str(J)])
disp(['Build Cost: ',num2str(Cost.C_build)])
disp(['Holding Cost: ',num2str(Cost.C_hold)])
disp(['Launch Cost: ',num2str(Cost.C_launch)])
disp(['Transfer Cost: ',num2str(Cost.C_transfer)])
disp(['P(Xi< N_sat): ',num2str(c(1)+p_loss_plane)])
disp(['P(Xp = 0): ',num2str(c(2)+p_loss_park)])