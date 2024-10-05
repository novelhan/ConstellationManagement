close all
clear all
clc

%% Test Param
%.. Sim time
dt_sim      =   45;                  % [day]
time_sim    =   0:dt_sim:365*100000;

%.. Marcov Chain Period
dt_period   =   90;                 % [day]
cnt_period  =   round(dt_period/dt_sim);

%.. Failure rate
p_fail      =   0.5/365;         % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_period    =   p_fail * dt_period;  % [#/dt_drfit]

%.. (Q,R) Policy Parameter
Q1  =   3;
R1  =   7;
Q2  =   2;
R2  =   5;

%.. Direct LV Parameters
mu_LV   =   40;     % [day]
T0_LV   =   0;
p_av    =   0.8;

%% (Q, R) Policy with fixed time order
rng('default')
iter_max    =   1000000;

Ni_Q2 = zeros(iter_max, 1);
idx_Q2 = zeros(iter_max, 1);
dT_Q2 = zeros(iter_max, 1);
cnt_Q2 = zeros(iter_max, 1);

xmax = Q1 + Q2 + R1;
xx = 0:1:xmax;
xx_edge = -0.5:1:(xmax+0.5);
cnt_edge = 0.5:1:(cnt_period+0.5);

for iter = 1:iter_max
    % Initialize Variables    
    N_on = 3;
    dt_LV = T0_LV + exprnd(mu_LV,1);
    cnt_LV = ceil(dt_LV/dt_sim);
    
    i = 0;
%     if rand < 0.1
%         i = 0;
%     else
%         i = 1;
%     end
    while 1
        k = mod(i,cnt_period) ;
        
        % 1. Check Q2 Resupply Arrival
        if cnt_LV == 1
            N_on = N_on + Q2;
            Ni_Q2(iter) = N_on;
            idx_Q2(iter) = k;
            dT_Q2(iter) = dt_LV;
            cnt_Q2(iter) = ceil(dt_LV/dt_sim);
            break
        elseif cnt_LV > 1
            cnt_LV = cnt_LV - 1;
        end
        
        % 2. Generate Failed Sample
        N_fail = poissrnd((Q1+R1)*p_sim, 1);
        N_on = max([N_on - N_fail, 0]);
        
        % 3. Check Q1 Resupply
        if k == 1
            if N_on <= R1 && rand < p_av
                N_on = N_on + Q1;
            end
        end
        
        i = i+1;
    end
end

%%

figure(1)
histogram(Ni_Q2(:),'Normalization','probability')

figure(2)
histogram(idx_Q2(:),'Normalization','probability')

figure(3)
histogram(dT_Q2(:),'Normalization','probability')


%%
f = p_sim;
lam = 1/mu_LV;

%.. State
xmax = Q1 + Q2 + R1;
x = xmax:-1:0;
nx = length(x);

%.. Failure Transition Matrix
Pf = zeros(nx);
for i = 1:nx
    Pf(i:end,i) = poisspdf(0:nx-i, (Q1+R1)*f)';        % Constant Failure Rate
    Pf(end,i) = 1 - sum(Pf(1:end-1, i));
end
i_r1 = xmax - R1;
i_r2 = xmax - R2;
U0 = Pf(1:i_r2, 1:i_r2);
U1 = Pf(i_r2+1:end, i_r2+1:end);
B = Pf(i_r2+1:end, 1:i_r2);

%.. Resupply Transition Matrix
Pq1 = zeros(nx);
Pq1(1:i_r1, 1:i_r1) = eye(Q1+Q2);
Pq1(i_r1+1-Q1:nx-Q1, i_r1+1:nx) = p_av*eye(R1+1);
Pq1(i_r1+1:nx, i_r1+1:nx) = Pq1(i_r1+1:nx, i_r1+1:nx) + (1-p_av)*eye(R1+1);

Pq2 = zeros(nx);
Pq2(1:Q2, 1:Q2) = eye(Q2);
Pq2(1:nx-Q2, Q2+1:nx) = eye(nx-Q2);

%.. Selection Matrix
C1 = eye(nx);
C2 = zeros(nx);
C2(i_r2+1:end,i_r2+1:end) = eye(R2+1);
C1 = C1 - C2;

Pr2o_11 = C2*Pq1*Pf*inv(eye(nx) - C1*Pf*C1*Pq1*Pf);
Pr2o_21 = C2*Pf*C1*Pq1*Pf*inv(eye(nx) - C1*Pf*C1*Pq1*Pf);
Pr2o_12 = C2*Pq1*Pf*C1*Pf*inv(eye(nx) - C1*Pq1*Pf*C1*Pf);
Pr2o_22 = C2*Pf*inv(eye(nx) - C1*Pq1*Pf*C1*Pf);

exp_dt = exp(-lam*dt_sim);
Po2r_11 = exp_dt*(1 - exp_dt)*Pq2*Pf*inv(eye(nx) - exp_dt^2*Pq1*Pf*Pf);
Po2r_21 = (1 - exp_dt)*Pq2*inv(eye(nx) - exp_dt^2*Pq1*Pf*Pf);
Po2r_12 = (1 - exp_dt)*Pq2*inv(eye(nx) - exp_dt^2*Pf*Pq1*Pf);
Po2r_22 = exp_dt*(1 - exp_dt)*Pq2*Pq1*Pf*inv(eye(nx) - exp_dt^2*Pf*Pq1*Pf);

Pr2o = [Pr2o_11, Pr2o_12;
        Pr2o_21, Pr2o_22];
Po2r = [Po2r_11, Po2r_12;
        Po2r_21, Po2r_22];
Po2o = Po2r*Pr2o;

xi = zeros(nx,1);
xi(10) = 1;
a = 1.0*Po2r_11*xi + 0.0*Po2r_12*xi;
b = 1.0*Po2r_21*xi + 0.0*Po2r_22*xi;

figure(1); hold on
plot(x,a+b,'r*')

figure(4); hold on
histogram(Ni_Q2((2*dt_sim < dT_Q2) & (dT_Q2 < 3*dt_sim)),'Normalization','probability')
plot(x,Pq2*Pq1*Pf*Pf*xi,'r*')


figure(5); hold on
histogram(cnt_Q2(:),'Normalization','probability')


ii = 0:10;
rho = (1 - exp_dt)*exp(-lam*dt_sim*ii)