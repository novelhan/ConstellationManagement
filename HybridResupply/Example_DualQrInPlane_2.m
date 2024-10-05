close all
clear all
clc
% T 구간에서 T/2 를 2번가는 것과, T를 한번가는것 그리고 T/n을 n번 가는것 같은지 확인
% 포아송, 무한스테이트면 같을 것 같은데 한정스테이트, 스테이트 의존 고장모델도 같은지 확인 아마 아닐듯
%% Test Param
%.. Sim time
dt_sim      =   45;                  % [day]
time_sim    =   0:dt_sim:365*5000;

%.. Marcov Chain Period
dt_period   =   90;                 % [day]
cnt_period  =   round(dt_period/dt_sim);

%.. Failure rate
p_fail      =   2.5/365;            % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_period    =   p_fail * dt_period; % [#/dt_drfit]

%.. (Q,R) Policy Parameter
Q1  =   3;
R1  =   4;
Q2  =   4;
R2  =   3;

%.. Direct LV Parameters
mu_LV   =   40;     % [day]
T0_LV   =   0;
p_av    =   0.8;

% %.. Failure rate
% p_fail      =   0.2/365;         % [#/day]
% p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
% p_period    =   p_fail * dt_period;  % [#/dt_drfit]
% 
% %.. (Q,R) Policy Parameter
% Q1  =   4;
% R1  =   40;
% Q2  =   6;
% R2  =   30;
% 
% %.. Direct LV Parameters
% mu_LV   =   40;     % [day]
% T0_LV   =   0;
% p_av    =   0.8;

%% (Q, R) Policy with fixed time order
rng('default')
iter_max    =   10;

xmax = Q1 + Q2 + R1;
xx = 0:1:xmax;
xx_edge = -0.5:1:(xmax+0.5);
cnt_edge = 0.5:1:(cnt_period+0.5);

for iter = 1:1
    % Initialize Variables
    N_on =   R1+round(Q1*rand);
    Ni_on =   zeros(length(time_sim),1);    % The number of availalbe stock at every moment
%     Ni_Q1 = -ones(length(time_sim),1);      % The number of stock right after the Q1 resupply moment
%     Ni_R1 = -ones(length(time_sim),1);      % The number of stock at the R1 reordering moment
%     Ni_Q2 = -ones(length(time_sim),1);      % The number of stock right after the Q2 resupply moment
%     Ni_R2 = -ones(length(time_sim),1);      % The number of stock at the R2 reordering moment
    Ni_Q1 = zeros(xmax+1,cnt_period);
    Ni_R1 = zeros(xmax+1,cnt_period);
    Ni_Q2 = zeros(xmax+1,cnt_period);
    Ni_R2 = zeros(xmax+1,cnt_period);
    T_LV = -ones(length(time_sim),1);       % Lead Time
    idx_Q2 = zeros(1,cnt_period);
    idx_R2 = zeros(1,cnt_period);
    N_av = 0;
    N_nav = 0;
    cnt_LV = 0;
    
    % Apply Policy
    for i = 1:length(time_sim)
        k = mod(i,cnt_period) + 1;
        
        % 1. Check Q2 Resupply Arrival
        if cnt_LV == 1
            N_on = N_on + Q2;
            cnt_LV = 0;
            Ni_Q2(N_on+1,k) = Ni_Q2(N_on+1,k) + 1;
            idx_Q2(1,k) = idx_Q2(1,k) + 1;
        elseif cnt_LV > 1
            cnt_LV = cnt_LV - 1;
        end
        
        % 2. Generate Failed Sample
%         N_fail = poissrnd(N_on*p_period, 1);
        N_fail = poissrnd((Q1+R1)*p_sim, 1);
        N_on = max([N_on - N_fail, 0]);
        
        % 3. Check Q1 Resupply
        if k == 1
            Ni_R1(N_on+1,k) = Ni_R1(N_on+1,k) + 1;
            if N_on <= R1
                if rand < p_av
                    N_on = N_on + Q1;
                    N_av = N_av + 1;
                else
                    N_nav = N_nav + 1;
                end
            end
            Ni_Q1(N_on+1,k) = Ni_Q1(N_on+1,k) + 1;
        end
        
        % 4. Check Q2 Resupply
        if cnt_LV == 0 && N_on <= R2
            dt_LV = T0_LV + exprnd(mu_LV,1);
            cnt_LV = ceil(dt_LV/dt_sim);
            T_LV(i) = cnt_LV*dt_sim;
            Ni_R2(N_on+1,k) = Ni_R2(N_on+1,k) + 1;
            idx_R2(1,k) = idx_R2(1,k) + 1;
        end
        
        % Update Profile
        Ni_on(i) = N_on;
    end
end
T_LV = T_LV(T_LV >= 0);

%%
figure(1)
plot(time_sim, Ni_on)

figure(2); hold on
histogram(Ni_on(:),'Normalization','probability')
xlabel('Number of stock for entire period')
ylabel('Probability')

figure(3); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Ni_Q1,2), 'Normalization','probability')
xlabel('Number of stock right after Q1 resupply')
ylabel('Probability')

figure(4); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Ni_R1,2), 'Normalization','probability')
xlabel('Number of stock at reordering of R1 resupply')
ylabel('Probability')

figure(5); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Ni_Q2,2), 'Normalization','probability')
xlabel('Number of stock right after Q2 resupply')
ylabel('Probability')

figure(6); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Ni_R2,2), 'Normalization','probability')
xlabel('Number of stock at reordering of R2 resupply')
ylabel('Probability')

% figure(7)
% histogram(T_LV(:),'Normalization','probability','BinWidth',1)
% xlabel('Lead Time Distribution')
% ylabel('Probability')
% 
% figure(51); hold on; grid on
% histogram('BinEdges', cnt_edge, 'BinCounts', idx_Q2, 'Normalization','probability')
% xlabel('Number of parking stock during entire period')
% ylabel('Probability')
% 
% figure(61); hold on; grid on
% histogram('BinEdges', cnt_edge, 'BinCounts', idx_R2, 'Normalization','probability')
% xlabel('Number of parking stock during entire period')
% ylabel('Probability')


mu_Q2 = idx_Q2/sum(idx_Q2);
mu_R2 = idx_R2/sum(idx_R2);
pi_Q2 = zeros(size(Ni_Q2));
pi_R2 = zeros(size(Ni_R2));
for i = 1:cnt_period
    pi_Q2(:,i) = Ni_Q2(:,i)/sum(Ni_Q2(:,i));
    pi_R2(:,i) = Ni_R2(:,i)/sum(Ni_R2(:,i));
end

% figure(52)
% hold on; grid on
% bar(xx, pi_Q2)
% 
% figure(62)
% hold on; grid on
% bar(xx, pi_R2)
% 
% figure(5)
% bar(xx, pi_Q2*mu_Q2')
% 
% figure(6)
% bar(xx, pi_R2*mu_R2')

%%

f = p_sim;
% f = p_period;
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
Po2o = Po2r*Pr2o;
Pr2r = Pr2o*Po2r;
% pi_o = limitdist(Po2o')         % Prob. Dist after reorder arrives
% pi_r = Pr2o*pi_o                % Prob. Dist when reorder is made
% pi_r = limitdist(Pr2r')


pi_Q2 = flip(pi_Q2);
pi_R2 = flip(pi_R2);
for i = 1:cnt_period
    pi_Q2(:,i) = mu_Q2(i)*pi_Q2(:,i);
    pi_R2(:,i) = mu_R2(i)*pi_R2(:,i);
end

%%
 
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
pi_o = limitdist(Po2o');
pi_Q2a = reshape(pi_o,nx,cnt_period);

P1_0 = (C1*Pq1*Pf + C1*Pf*C1*Pq1*Pf)*inv(eye(nx) - C1*Pf*C1*Pq1*Pf);
P1_1 = (C1*Pf + C1*Pq1*Pf*C1*Pf)*inv(eye(nx) - C1*Pq1*Pf*C1*Pf);
pi_1 = P1_0*pi_Q2a(:,1) + P1_1*pi_Q2a(:,2);

pi_r = Pr2o*pi_o;
pi_R2a = reshape(pi_r,nx,cnt_period);

P2_0 = (eye(nx) + exp_dt*Pf)*inv(eye(nx) - exp_dt^2*Pq1*Pf*Pf);
P2_1 = (eye(nx) + exp_dt*Pq1*Pf)*inv(eye(nx) - exp_dt^2*Pf*Pq1*Pf);
pi_2 = P2_0*pi_R2a(:,1) + P2_1*pi_R2a(:,2);

pi_t = 1/(sum(pi_1)+sum(pi_2))*pi_1 + 1/(sum(pi_1)+sum(pi_2))*pi_2;

% figure(2)
% plot(x, pi_t, 'r*')
% 
% figure(5)
% plot(x, sum(pi_Q2a,2), 'r*')
% 
% figure(6)
% plot(x, sum(pi_R2a,2), 'r*')

rho = (1 - exp_dt);
A11 = Pf*C1 + rho*Pf*Pq2*C2;
A12 = rho*Pf*Pq2;
A21 = (1-rho)*Pf*C2;
A22 = (1-rho)*Pf;
AA = [A11, A12; A21, A22];
P_AA = AA^2 * blkdiag(Pq1,Pq1);
pi_AA = limitdist(P_AA');

pi_f0 = pi_AA;
pi_f1 = AA * blkdiag(Pq1,Pq1)*pi_f0;
pi_r0 = C2*Pq1*pi_f0(1:12)/sum(C2*Pq1*pi_f0(1:12));
pi_r1 = C2*pi_f1(1:12)/sum(C2*pi_f1(1:12));
pi_R2A = (C2*Pq1*pi_f0(1:12) + C2*pi_f1(1:12))/sum(C2*Pq1*pi_f0(1:12) + C2*pi_f1(1:12));
pi_tt = (Pq1*pi_f0(1:12) + pi_f1(1:12) + Pq1*pi_f0(13:end) + pi_f1(13:end))/2;
pi_q1_before = pi_f0(1:12) + pi_f0(13:end);
pi_q1_after = Pq1*pi_q1_before;

figure(2)
plot(x, pi_tt, 'r*')

figure(3)
plot(x, pi_q1_after, 'r*')

figure(4)
plot(x, pi_q1_before, 'r*')

% figure(5)
% plot(x, sum(pi_Q2a,2), 'r*')

figure(6)
plot(x, pi_R2A, 'r*')

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