clear all
close all
clc

%% System Parameter
%.. Marcov Chain Period
dt_mc   =   1;                 % [day]
PARA.dt_mc = dt_mc;

%.. Failure rate
% 0.05, 0.1, 0.2
f_fail = 0.2/365; % [#/day]
f_type = 1; % 0:Const, 1:State-Dependant
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
mu_LV   =   60;     % [day]
T0_LV   =   30;
PARA.mu_LV = 1/mu_LV;
PARA.T0_LV = T0_LV;

%.. Desired System Performance
n_sat   =   40;
p_loss  =   0.05;
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
opts = optimoptions('ga','MaxStallGenerations',50,'FunctionTolerance',1e-10,...
                    'MaxGenerations',500,'PopulationSize',400,'Display','iter');
% opts = optimoptions('Display','iter','PopulationSize',400);
rng default % For reproducibility
% x_opt = ga(@(x) CostFun(x,PARA), 2, [], [], [], [], [Qmin Rmin], [Qmax Rmax],...
%            @(x) ConstFun(x,PARA), intcon, opts);

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
    J_hold_i = N_plane * ParaCost.c_holding * sum(xx_i(1:idx_i).*PI_i.pi_ir(1:idx_i));
    J_hold_p = N_park * ParaCost.c_holding * Qi * sum(xx_p(1:idx_p).*PI_p.pi_dr(1:idx_p));
    
    %.. Total Cost
    J = J_build + J_launch + J_hold_i + J_hold_p;
end

%% Constraint Function
function [c,ceq] = ConstFun(x, ParaCost, ParaInPlane, ParaParking)
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

    %.. 

    %.. Qp >= Rp + 1
    c(2) = ParaParking.Q - ParaParking.R - 1;
    
    %.. Qi*Qp <= Q_max TODO
    % c(3) = x(1) * x(3) - ParaCost.Q_max;
    
    %.. 0 < Rp < Qp_max (Handled by external feature)
    %.. 

    
    idx = find(xx == PARA.n_sat);
    p_loss = sum(PI.pi_dr(idx+1:end));
    c = p_loss - PARA.p_loss;
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