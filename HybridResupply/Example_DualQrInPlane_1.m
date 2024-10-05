close all
clear all
clc

% This script is for testing dual source replenishment strategy
% with assumption of synchronized resupply time

%% Test Param
%.. Sim time
dt_sim      =   90;                  % [day]
time_sim    =   0:dt_sim:365*100000;

%.. Marcov Chain Period
dt_period   =   90;                 % [day]
cnt_period  =   round(dt_period/dt_sim);

%.. Failure rate
p_fail      =   0.05/365;            % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_period    =   p_fail * dt_period; % [#/dt_drfit]

%.. (Q,R) Policy Parameter
Q1  =   4;
R1  =   42;
Q2  =   4;
R2  =   38;

%.. Direct LV Parameters
mu_LV   =   60;     % [day]
T0_LV   =   0*dt_period;
p_av    =   0.9;

%% (Q, R) Policy with fixed time order
rng('default')
iter_max    =   10;

for iter = 1:1
    % Initialize Variables
    N_on =   R1+round(Q1*rand);
    Ni_on =   zeros(length(time_sim),1);    % The number of availalbe stock at every moment
    Ni_Q1 = -ones(length(time_sim),1);      % The number of stock right after the Q1 resupply moment
    Ni_R1 = -ones(length(time_sim),1);      % The number of stock at the R1 reordering moment
    Ni_Q2 = -ones(length(time_sim),1);      % The number of stock right after the Q2 resupply moment
    Ni_R2 = -ones(length(time_sim),1);      % The number of stock at the R2 reordering moment
    T_LV = -ones(length(time_sim),1);       % Lead Time
    N_av = 0;
    N_nav = 0;
    Lv_on = 0;
    cnt = 0;
    dcnt = 0;
    % Apply Policy
    for i = 1:length(time_sim)

        % Check Resupply Arrival
        if Lv_on == 1
            N_on = N_on + Q2;
            Lv_on = 0;
            Ni_Q2(i) = N_on;
        elseif Lv_on > 1
            Lv_on = Lv_on - 1;
        end
        
        % Generate Fail Sample at every drift period
        if mod(i,cnt_period) == 0
%             N_fail = poissrnd(N_on*p_period, 1);
            N_fail = poissrnd((Q1+R1)*p_period, 1);
            N_on = max([N_on - N_fail, 0]);
        end
        
        % Check Reorder at every drift period
        if mod(i,cnt_period) == 0
            % Check Q1 First
            
            Ni_R1(i) = N_on;
            if N_on <= R1
                if rand < p_av
                    N_on = N_on + Q1;
                    N_av = N_av + 1;
                else
                    N_nav = N_nav + 1;
                end
            end
            Ni_Q1(i) = N_on;
            % Check Q2 Next
            if Lv_on == 0 && N_on <= R2
                dt_LV = T0_LV + exprnd(mu_LV,1);
                Lv_on = ceil(dt_LV/dt_sim);
                T_LV(i) = Lv_on - 1;
                Ni_R2(i) = N_on;
            end
            
            if Lv_on ~= 0
                cnt = cnt + 1;
            else
                dcnt = dcnt + 1;
            end
        end
        
        % Update Profile
        Ni_on(i) = N_on;
    end
end
% Ni_on = Ni_on(Ni_on >= 0);
Ni_Q1 = Ni_Q1(Ni_Q1 >= 0);
Ni_R1 = Ni_R1(Ni_R1 >= 0);
Ni_Q2 = Ni_Q2(Ni_Q2 >= 0);
Ni_R2 = Ni_R2(Ni_R2 >= 0);
T_LV = T_LV(T_LV >= 0);

%%
figure(1)
plot(time_sim, Ni_on)

figure(2)
histogram(Ni_on(cnt_period:cnt_period:end),'Normalization','probability')
xlabel('Number of stock for entire period')
ylabel('Probability')

figure(3)
histogram(Ni_Q1(:),'Normalization','probability')
xlabel('Number of stock right after Q1 resupply')
ylabel('Probability')

figure(4)
histogram(Ni_R1(:),'Normalization','probability')
xlabel('Number of stock at reordering of R1 resupply')
ylabel('Probability')

figure(31)
histogram(Ni_Q2(:),'Normalization','probability')
xlabel('Number of stock right after Q2 resupply')
ylabel('Probability')

figure(41)
histogram(Ni_R2(:),'Normalization','probability')
xlabel('Number of stock at reordering of R2 resupply')
ylabel('Probability')

figure(5)
histogram(T_LV(:),'Normalization','probability','BinWidth',1)
xlabel('Lead Time Distribution')
ylabel('Probability')

%%
% f = p_sim;
f = p_period;
lam = 1/mu_LV;

%.. State
xmax = Q1 + Q2 + R1;
x = xmax:-1:0;
nx = length(x);

%.. Failure Transition Matrix
Pf = zeros(nx);
for i = 1:nx
    Pf(i:end,i) = poisspdf(0:nx-i, (Q1+R1)*f)';        % Constant Failure Rate
%     Pf(i:end,i) = poisspdf(0:nx-i, (nx-i)*f)';  % Stock Level Dependant Failure Rate
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

%.. Full Transition
% Pr2o = C2*inv(eye(nx) - Pq1*Pf*C1);
Pr2o = C2*Pq1*Pf*inv(eye(nx) - C1*Pq1*Pf);
% Po2r = (1 - exp(-lam*dt_period))*Pq2*Pq1*Pf*inv(eye(nx)-exp(-lam*dt_period)*Pq1*Pf); 
Po2r = (1 - exp(-lam*dt_period))*Pq2*inv(eye(nx)-exp(-lam*dt_period)*Pq1*Pf);
Pr2r = Po2r*Pr2o;
Po2o = Pr2o*Po2r;
pi_r = limitdist(Pr2r');         % Prob. Dist after reorder arrives
pi_o = Pr2o*pi_r;                % Prob. Dist when reorder is made
% pi_r = limitdist(Pr2r')

pi_1f = Pf*inv(eye(nx) - C1*Pq1*Pf)*pi_r;
pi_2f = exp(-lam*dt_period)*Pf*inv(eye(nx) - exp(-lam*dt_period)*Pq1*Pf)*pi_o;
pi_r1 = 1/(sum(pi_1f)+sum(pi_2f))*pi_1f + 1/(sum(pi_1f)+sum(pi_2f))*pi_2f;
% pi_r1 = 0.5/(sum(pi_1f))*pi_1f + 0.5/(sum(pi_2f))*pi_2f;
pi_q1 = Pq1*pi_r1;

pi_1 = C1*Pq1*Pf*inv(eye(nx) - C1*Pq1*Pf)*pi_r;
pi_2 = inv(eye(nx) - exp(-lam*dt_period)*Pq1*Pf)*pi_o;
pi_t = 1/(sum(pi_1)+sum(pi_2))*pi_1 + 1/(sum(pi_1)+sum(pi_2))*pi_2;

figure(2); hold on
plot(x, pi_t, 'r*')

figure(3); hold on
plot(x, pi_q1, 'r*')
figure(4); hold on
plot(x, pi_r1, 'r*')

figure(31); hold on
plot(x, pi_r, 'r*')

figure(41); hold on
plot(x, pi_o, 'r*')

rho1 = 1 - exp(-lam*dt_period);
a = Pf*C1*Pq1 + rho1*Pf*Pq2*inv(eye(nx) - exp(-lam*dt_period)*Pq1*Pf)*C2*Pq1;
pi_0 = limitdist(a');

% pi_1 = inv(eye(nx) - C1*Pq1*Pf)*C1*Pq1*pi_0
% pi_1 = C1*Pq1*inv(eye(nx) - C1*Pq1*Pf)*pi_0
% pi_2 = inv(eye(nx) - exp(-lam*dt_period)*Pq1*Pf)*C2*Pq1*pi_0
% pi_2 = inv(eye(nx) - exp(-lam*dt_period)*Pq1*Pf)*C2*Pq1*pi_0

%% Frequency based method
%.. No Wait Time
rho = (1 - exp(-lam*dt_period));
A11 = Pf*C1 + rho*Pf*Pq2*C2;
A12 = rho*Pf*Pq2;
A21 = (1-rho)*Pf*C2;
A22 = (1-rho)*Pf;
AA = [A11, A12; A21, A22];
P_AA = AA * blkdiag(Pq1,Pq1);
pi_AA = limitdist(P_AA');

%.. Const Wait Time + Exp. Dist
% rho = (1 - exp(-lam*dt_period));
% A11 = Pf*C1*Pq1;
% A12 = rho*Pf*Pq2*Pq1;
% A21 = Pf*Pq1*Pf*C2*Pq1;
% A22 = (1-rho)*Pf*Pq1;
% AA = [A11, A12; A21, A22];
% pi_AA = limitdist(AA');
pi_fo = C2*Pq1*pi_AA(1:nx)/sum(C2*Pq1*pi_AA(1:nx));
pi_ft = Pq1*sum(reshape(pi_AA,nx,2),2);
tmp1 = reshape(pi_AA,nx,2);
tmp1 = [tmp1, Pf*C2*Pq1*tmp1(:,1)];

% pi_ft = Pq1*sum(tmp1,2)/sum(sum(tmp1,2));

% rho = (1 - exp(-lam*dt_period));
% A11 = Pf*C1*Pq1;
% A12 = rho*Pf*Pq2*Pq1;
% A13 = zeros(nx,nx);
% A21 = zeros(nx,nx);
% A22 = (1-rho)*Pf*Pq1;
% A23 = Pf*Pq1;
% A31 = Pf*C2*Pq1;
% A32 = zeros(nx,nx);
% A33 = zeros(nx,nx);
% AA = [A11, A12, A13; 
%       A21, A22, A23;
%       A31, A32, A33];
% pi_AA = limitdist(AA');
% pi_fo = C2*Pq1*pi_AA(1:nx)/sum(C2*Pq1*pi_AA(1:nx));
% pi_ft = Pq1*sum(reshape(pi_AA,nx,3),2);
% tmp2 = reshape(pi_AA,nx,3)

figure(2)
plot(x, pi_ft, 'gd')

figure(41)
plot(x, pi_fo, 'gd')

function p=limitdist(P)
%Obtain the stationary probability distribution
%vector p of an irreducible, recurrent Markov
%chain by state reduction. P is the transition
%probabilities matrix of a discrete-time Markov
%chain or the generator matrix Q.
% https://www.math.wustl.edu/~feres/Math450Lect04.pdf

[ns ms]=size(P);
n=ns;
p=zeros(n);
while n>1
    n1=n-1;
    s=sum(P(n,1:n1));
    P(1:n1,n)=P(1:n1,n)/s;
    n2=n1;
    while n2>0
        P(1:n1,n2)=P(1:n1,n2)+P(1:n1,n)*P(n,n2);
        n2=n2-1;
    end
    n=n-1;
end
%backtracking
p(1)=1;
j=2;
while j<=ns
    j1=j-1;
    p(j)=sum(p(1:j1).*(P(1:j1,j))');
    j=j+1;
end
p=p/(sum(p));
end