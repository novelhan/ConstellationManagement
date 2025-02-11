% This script optimizes the direct spare management strategy using the genetic algorithm.
% The decision variables to optimize are (r,q) of the inventory control
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

%% System Parameter
%.. Marcov Chain Period
dt_mc   =   1;                 % [day]

%.. Failure rate
f_ref = 0.1/365; % [#/day]
f_type = 1; % 0:Constant failure, 1:State Dependent failure

%.. Analysis Method
method = 0; % 0:Time based, 1:Ratio based

%.. Constellation Orbit Parameters (Given)
N_plane     =   40;                     % The number of orbital planes for constellation
N_sat       =   40;                     % The number of desigend satellites for each plane

%.. (small) Direct LV Parameters (goes to in-plane orbit)
mu_lv_small     =   30;
dt_lv_small     =   15;                % [day]

%.. Cost Model
c_build = 0.5; % Manufacturing cost per satellite [M$/sat]
c_hold = 0.5/365; % Holding cost per satellite in In-Plane orbit [M$/sat/day]
c_lv_small_part = 3; %[M$/sat/launch] for direct LV (small size)
c_lv_small_full = 7.5;%[M$/launch] for direct LV (small size)

%.. Desired System Performance
p_loss  =   0.05; % Resilience Prob, the prob. having X < n_sat must be smaller than p_loss

%.. Save the input structure
ParaInPlane.dt_mc = dt_mc;
ParaInPlane.f_ref = f_ref;
ParaInPlane.f_type = f_type;
ParaInPlane.N_plane = N_plane;
ParaInPlane.N_sat = N_sat;
ParaInPlane.mu_lv = mu_lv_small;
ParaInPlane.dt_lv = dt_lv_small;
ParaInPlane.method = method;

ParaCost.c_build = c_build;
ParaCost.c_hold = c_hold;
ParaCost.c_lv_small_part = c_lv_small_part;
ParaCost.c_lv_small_full = c_lv_small_full;

ParaConst.p_loss = p_loss;

%% Design Variable
%.. Launch Vehicle Size Model (Design Var.)
Qmin = 1;
Qmax = 3;

%.. Reorder Point Range (Design Var.)
Rmin = N_sat - 2*Qmax; % Considered to reduce the search space
Rmax = N_sat + 2*Qmax; 

%% Optimize
%.. Lower/Upper Bound
LB = [Qmin, Rmin];
UB = [Qmax, Rmax];

%.. Integer Index and Opt Option
intcon = [1 2];
opts = optimoptions('ga','MaxStallGenerations',20,'FunctionTolerance',1e-10,...
                    'MaxGenerations',50,'PopulationSize',100,'Display','iter');

%.. Run Opt
rng default % For reproducibility
x_opt = ga(@(x) CostDirectResupply(x,ParaCost,ParaInPlane), 2, [], [], [], [], LB, UB,...
           @(x) ConstDirectResupply(x,ParaConst,ParaInPlane), intcon, opts);
CostDirectResupply(x_opt,ParaCost,ParaInPlane)
%% Analysis Plot
qq = Qmin:Qmax;
rr = Rmin:Rmax;
[QQ,RR] = meshgrid(qq,rr);
Jset = zeros(length(qq),length(rr));
Cset = zeros(length(qq),length(rr));
feas = [];
infeas = [];
for i = 1:length(qq)
for j = 1:length(rr)
    x = [qq(i), rr(j)];
    Jset(i,j) = CostDirectResupply(x,ParaCost,ParaInPlane);
    Cset(i,j) = ConstDirectResupply(x,ParaConst,ParaInPlane);

    if Cset(i,j) >= 0
        infeas = [infeas; x, Jset(i,j)];
    else
        feas = [feas; x, Jset(i,j)];
    end
end
end
    
figure(1); hold on; grid on
surf(QQ',RR',Jset)
xlabel('Q')
ylabel('R')
zlabel('Cost')
% clim([min(feas(:,3)), max(feas(:,3))])
colorbar
[~,idx_min] = min(feas(:,3));
opt = feas(idx_min,:);
feas(idx_min,:) = [];
h1 = plot3(infeas(:,1), infeas(:,2), infeas(:,3)+0.1, 'rx', 'MarkerSize', 8);
h2 = plot3(feas(:,1), feas(:,2), feas(:,3)+0.1, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
h3 = plot3(opt(1), opt(2), opt(3)+0.1, 'y*', 'MarkerSize', 8);
legend([h1 h2 h3], 'Infeasible', 'Feasible', 'Optimal')

