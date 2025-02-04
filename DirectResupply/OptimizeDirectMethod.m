% This script optimizes the direct spare management strategy using the genetic algorithm.
% The decision variables to optimize are (r,q) of the inventory control
% The cost function is the effective operational cost of the strategy per unit cycle
% The constraints include the payload capacity of LV, resilience of the spare strategy
% The script requires optimization toolbox.

clear all
close all
clc

%% System Parameter
%.. Marcov Chain Period
dt_mc   =   1;                 % [day]
PARA.dt_mc = dt_mc;

%.. Failure rate
f_fail = 0.2/365; % [#/day] 0.05, 0.1, 0.15
f_type = 1; % 0:Constant failure, 1:State Dependent failure
PARA.f = f_fail;
PARA.f_type = 0;

%.. Launch Vehicle Cost Model
c_build = 0.5; %[$/sat]
c_holding = 0.5/365; %[$/sat/day]
c_lv_part = 10; %[$/sat/launch]
lv_full_discount = 0.05; %[%/launch]
PARA.c_build = c_build;
PARA.c_holding = c_holding;
PARA.c_lv_part = c_lv_part;
PARA.lv_full_discount = lv_full_discount;

%.. Launch Vehicle Lead Time Parameters
mu_LV   =   60; % Mean lead time [day]
T0_LV   =   30; % Constant shift lead time [day]
PARA.mu_LV = 1/mu_LV;
PARA.T0_LV = T0_LV;

%.. Desired System Performance
n_sat   =   40; % Number of operating satellite per orbit for nominal condition
p_loss  =   0.05; % Resilience Prob, the prob. having X < n_sat must be smaller than p_loss
PARA.n_sat = n_sat;
PARA.p_loss = p_loss;

%.. Launch Vehicle Size Model (Design Var.)
Qmin = 1;
Qmax = 6;
PARA.Qmax = Qmax;

%.. Reorder Point Range (Design Var.)
Rmin = n_sat;
Rmax = Rmin + 2*Qmax;


%% Optimize
x0 = [Qmax, Rmax];
intcon = [1 2];
opts = optimoptions('ga','MaxStallGenerations',20,'FunctionTolerance',1e-10,...
                    'MaxGenerations',50,'PopulationSize',100,'Display','iter');
rng default % For reproducibility
x_opt = ga(@(x) CostFun(x,PARA), 2, [], [], [], [], [Qmin Rmin], [Qmax Rmax],...
           @(x) ConstFun(x,PARA), intcon, opts);

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
    Jset(i,j) = CostFun(x,PARA);
    Cset(i,j) = ConstFun(x,PARA);

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
%%
function J = CostFun(x,PARA)
    Qdr = x(1);
    Rdr = x(2);
    [xx, PI, T] = ExactDirectProb(PARA.f, PARA.f_type, Qdr,Rdr,...
                                  PARA.mu_LV, PARA.dt_mc, PARA.T0_LV, PARA.n_sat);
    
    %.. Cost for making spare satellites per unit time
    J_build = PARA.c_build * Qdr / T.T_dr;
    
    %.. Cost for LV launching
    if Qdr == PARA.Qmax
        J_launch = (1 - PARA.lv_full_discount) * PARA.c_lv_part * PARA.Qmax / T.T_dr;
    else
        J_launch = PARA.c_lv_part * Qdr / T.T_dr;
    end
    
    %.. Cost for holding redundant spares
    idx = find(xx == PARA.n_sat + 1);
    J_hold = PARA.c_holding * sum(xx(1:idx).*PI.pi_dr(1:idx));
    
    %.. Total Cost
    J = J_build + J_launch + J_hold;
end

function [c,ceq] = ConstFun(x,PARA)
    Qdr = x(1);
    Rdr = x(2);
    [xx, PI, ~] = ExactDirectProb(PARA.f, PARA.f_type, Qdr,Rdr,...
                                  PARA.mu_LV, PARA.dt_mc, PARA.T0_LV, PARA.n_sat);
    
    idx = find(xx == PARA.n_sat);
    p_loss = sum(PI.pi_dr(idx+1:end));
    c = p_loss - PARA.p_loss;
    ceq = [];
end