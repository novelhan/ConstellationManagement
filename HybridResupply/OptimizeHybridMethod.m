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
ParaInPlane.dt_mc = dt_mc;
ParaParking.dt_mc = dt_mc;

%.. Failure rate
% 0.05, 0.1, 0.2
f_fail = 0.1/365; % [#/day]
f_type = 1; % 0:Const, 1:State-Dependant

%.. Constellation Orbit Parameters (Given)
N_plane     =   40;                     % The number of orbital planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane
h_plane     =   1200;                   % Altitude of in-plane orbits [km]
i_plane     =   50 * D2R;               % Inclination of in-plane orbits [rad]

%.. Parking Orbit Parameters (Given)
i_park      =   50 * D2R;               % Inclination of parking orbits [rad]

%.. Indirect LV Parameters (goes to parking orbit)
mu_lv_park      =   60;
dt_lv_park      =   30;                % [day]

%.. Direct LV Parameters (goes to in-plane orbit)
mu_lv_plane     =   30;
dt_lv_plane     =   15;                % [day]

%.. Cost Model
c_build = 0.5; % Manufacturing cost per satellite [M$/sat]
c_holding_plane = 0.5/365; % Holding cost per satellite in In-Plane orbit [M$/sat/day]
c_holding_park = 0.5/365; %Holding cost per satellite in Parking orbit [M$/sat/day] [M$/sat/day]
c_transfer = 0.5; % Risk and additional cost for indirect transfer [M$/batch]
c_fuel = 0.001; % Fuel cost for indirect transfer [M$/batch]
c_lv_plane_part = 3; %[M$/sat/launch] for direct LV
c_lv_plnae_full = 7.5;%[M$/launch] for direct LV
c_lv_park_part = 2; %[M$/sat/launch] for indirect LV
c_lv_park_full = 67; %[M$/launch] for indirect LV
lv_plane_max_capa = 3; %[#sat] max capacitiy for direct LV
lv_park_max_capa = 40; %[#sat] max capacitiy for indirect LV

ParaCost.c_build = c_build;
ParaCost.c_holding = c_holding_plane;
ParaCost.c_lv_part = c_lv_park_part;
ParaCost.lv_full_discount = lv_full_park_discount;
ParaCost.Qmax = Qp_park_max;

%.. Desired System Performance
p_loss  =   0.05; % Prob. having less than N_sat number of satellite
ParaConst.p_loss = p_loss;


%% Design Variable
%.. In-Plane InDirect Reorder Size Range (Design Var.)
Qi1_min = 1;
Qi1_max = 10;

%.. In-Plane Direct Reorder Size Range (Design Var.)
Qi2_min = 1;
Qi2_max = 4;



%.. In-Plane Reorder Point Range (Design Var.)
Ri_min = N_sat - 5;
Ri_max = N_sat + 5;

%.. Launch Vehicle Size Range (Design Var.)
Qp_min = 2;
Qp_park_max = Qp_park_max;

%.. Reorder Point Range (Design Var.)
Rp_min = 0;
Rp_max = Qp_park_max - 1;

%.. # of Parking Orbit Range
Np_min = 2;
Np_max = 20;

%.. Parking Orbit Altitude Range
hp_min = 500;
hp_max = 1000;


%% Optimize
%.. Lower/Upper Bound
LB = [Qi1_min, Ri_min, Qp_min, Rp_min, Np_min, hp_min];
UB = [Qi1_max, Ri_max, Qp_park_max, Rp_max, Np_max, hp_max];

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
    
    %.. Cost per unit time step for each different list
    C_build = ComputeManufacturingCost(T_dr, T_ir, ParaCost, ParaInPlane, ParaParking);
    C_hold = ComputeHoldingCost(P_Xi, P_Xp, ParaCost, ParaInPlane, ParaParking);
    C_launch = ComputeLaunchCost(T_dr, T_ir, ParaCost, ParaInPlane, ParaParking);
    C_transfer = ComputeOrbitTransferCost(P_Di, T_plane, ParaCost, ParaInPlane, ParaParking);
    
    %.. Total Cost
    J = C_build + C_hold + C_launch + C_transfer;
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

%% Cost function
function [C_build] = ComputeManufacturingCost(T_dr, T_ir, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    c_build = ParaCost.c_build;
    N_plane = ParaInPlane.N_plane;
    N_park = ParaParking.N_park;
    qi1 = ParaInPlane.Q1;
    qi2 = ParaInPlane.Q2;
    qp = ParaParking.Q;
    
    % Manufacturing Cost per unit time step
    C_build = c_build*(N_plane/T_dr*qi2 + N_park/T_ir*qi1*qp);
end

function [C_hold] = ComputeHoldingCost(P_Xi, P_Xp, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    c_hold = ParaCost.c_build;
    N_plane = ParaInPlane.N_plane;
    N_park = ParaParking.N_park;
    N_sat = ParaInPlane.N_sat;
    qi1 = ParaInPlane.Q1;
    
    % Avg. # of Spares
    Si_k = (length(P_Xi)-1:-1:0) - N_sat; % (Xi_max, ... , 1, 0) - N_sat
    Si_k(Si_k <= 0) = 0;
    Si_avg = dot(Si_k,P_Xi);
    
    Xp_k = length(P_Xp)-1:-1:0; % Xp_max, ... , 1, 0
    Xp_avg = dot(Xp_k,P_Xp);
    
    % Holding Cost per unit time step
    C_hold = c_hold*(N_plane*Si_avg + N_park*qi1*Xp_avg);
end

function [C_launch] = ComputeLaunchCost(T_dr, T_ir, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    c_launch_small_full = ParaCost.c_launch_small_full;
    c_launch_small_part = ParaCost.c_launch_small_part;
    c_launch_heavy_full = ParaCost.c_launch_heavy_full;
    c_launch_heavy_part = ParaCost.c_launch_heavy_part;
    N_plane = ParaInPlane.N_plane;
    N_park = ParaParking.N_park;
    qi1 = ParaInPlane.Q1;
    qi2 = ParaInPlane.Q2;
    qp = ParaParking.Q;
    
    % Apply minimum cost among full / partial cost
    c_launch_small = min(c_launch_small_full, qi2*c_launch_small_part);
    c_launch_heavy = min(c_launch_heavy_full, qi1*qp*c_launch_heavy_part);
    
    % Launch Cost per unit time step
    C_launch = (N_plane/T_dr)*c_launch_small + (N_park/T_ir)*c_launch_heavy;
end


function [C_trn] = ComputeOrbitTransferCost(P_Di, T_plane, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    h_park = ParaParking.h;
    h_plane = ParaInPlane.h;
    N_plane = ParaInPlane.N_plane;
    m_sat = ParaInPlane.m_sat;
    m_bus = ParaParking.m_bus;
    Vex = ParaParking.Vex;
    q1i = ParaInPlane.Q1;
    c_fuel = ParaCost.c_fuel;
    c_trn = ParaCost.c_trn;
    
    % Fuel for both (bus + spares) transfer
    m_fuel = ComputeTransferFuelPerBatch(h_park, h_plane, m_bus, m_sat, q1i, Vex);
    
    % Avg. # of transfer for every review period (T_plane)
    Di_k = 0:length(P_Di)-1; % 0,1,...,Dmax
    Di_avg = dot(Di_k,P_Di);
    
    % Transfer cost per unit time step
    C_trn = N_plane/T_plane * Di_avg * ( c_fuel*m_fuel + c_trn );
end

function [m_fuel] = ComputeTransferFuelPerBatch(h_park, h_plane, m_bus, m_sat, q1i, Vex)
    global R_earth mu_earth
    ai = R_earth + h_plane;
    ap = R_earth + h_park;
    DV = sqrt(mu_earth/ap)*( sqrt(2*ai/(ai+ap)) - 1 ) + sqrt(mu_earth/ai)*( 1 - sqrt(2*ap/(ai+ap)) );
    m_fuel = (m_bus + q1i*m_sat)*( exp(DV/Vex) - 1 );
end


%% RAAN Period Computing Function
function [dt_plane, dt_park] = ComputeRaanPeriod(ParaInPlane, ParaParking)
    %.. Earth Parameters
    global R_earth mu_earth J2_earth
    
    %.. In-Orbit Parameters
    a_plane     =   R_earth + ParaInPlane.h;
    e_plane     =   0;
    i_plane     =   ParaInPlane.inc;
    n_plane     =   sqrt(mu_earth/a_plane^3);
    N_plane     =   ParaInPlane.N_orbit;
    Wdot_orbit  =   -3/2 * J2_earth* n_plane * R_earth^2 / (a_plane*(1-e_plane))^2 * cos(i_plane);
    
    %.. Parking-Orbit Parameters
    a_park      =   R_earth + ParaParking.h;
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