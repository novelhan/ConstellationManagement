close all
clear all
clc
% T 구간에서 T/2 를 2번가는 것과, T를 한번가는것 그리고 T/n을 n번 가는것 같은지 확인
% 포아송, 무한스테이트면 같을 것 같은데 한정스테이트, 스테이트 의존 고장모델도 같은지 확인 아마 아닐듯
%% Test Param
%.. Sim time
dt_sim      =   1;                  % [day]
time_sim    =   0:dt_sim:365*500;

%.. Marcov Chain Period
dt_mc       =   dt_sim;             % [day]

%.. In-plane Period (Time duration of in-plane for subsequent parking contact, RAAN Drift time)
dt_plane   =   90;                 % [day]
cnt_plane  =   round(dt_plane/dt_sim);

%.. Failure rate
p_fail      =   0.15/365;           % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_mc        =   p_fail * dt_mc;     % [#/dt_mc]
p_type      =   0; % 0 for const, 1 for state dependant

%.. Hybrid In-plane (Q1,R1,Q2,R2) Policy Parameter
Q1  =   3;
R1  =   8;
Q2  =   2;
R2  =   4;
n_sat   =   40; % Nominal Satellite Count

%.. State Parameter
Xmax = max(Q1+Q2+R1, Q1+Q2+R2); % Max State Level
Xnum = 0:1:Xmax; % State Counts

%.. Parking Availablity (Test distribution)
Dmax = ceil(Xmax/Q1);
Pav = sqrt(Dmax+1:-1:1);
Pav = Pav/sum(Pav);
Pav_sum = cumsum(Pav);
Kappa = [1, 1 - Pav_sum(1:end-1)];
Nav = 0:1:Dmax;

%.. Direct LV Parameters
mu_LV   =   40;     % [day]
T0_LV   =   0;
dTlv_max = ceil((T0_LV + mu_LV + 5*mu_LV)/dt_mc); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]
dTlv = 1:dTlv_max; % Lead Time Bin (Minimum is set as dt_sim)

%% (Q1, R1, Q2, R2) Policy
rng('default')
iter_max = 1; % Number of different initial condition

% The number of availalbe stock at each time step / histogram
Non = zeros(length(time_sim), iter_max);
Xon = zeros(Xmax+1, iter_max);

% The histogram of the number of stock for the indirect resupply moment
Xq1 = zeros(Xmax+1, iter_max);
Xr1 = zeros(Xmax+1, iter_max);

% The histogram of the number of stock for the direct resupply moment
Xq2 = zeros(Xmax+1, cnt_plane, iter_max);
Xr2 = zeros(Xmax+1, cnt_plane, iter_max);
Idxq2 = zeros(cnt_plane, iter_max); % Time index counter for q2 arrival
Idxr2 = zeros(cnt_plane, iter_max); % Time index counter for q2 reorder

% Direct Lead time distribution
Xlv = zeros(dTlv_max, iter_max); 

% Parking available histogram
Xav = zeros(Dmax+1, iter_max);  
Xdmd = zeros(Dmax+1, iter_max); 

%% Run Each Simulation
for iter = 1:iter_max
    % The number of satellite at current time step
    Non_k = Xmax - round(Q1*rand) - round(Q2*rand);
    
    % ETC 
    cnt_lv = -1; % -1 for not ordered
    lv_cnt = 0;
    
    % Apply Policy
    for k = 1:length(time_sim)
        % RAAN Contact Counter
        cnt_p = mod(k,cnt_plane) + 1;
        
        %%% 1. Generate Fail Sample at Every Contact
        if p_type == 0 %.. Const Failure Rate
            N_fail = CustomPoisRnd(n_sat*p_sim, 1);
        else %.. State Dependant Failure Rate
            if Non_k > n_sat
                N_fail = CustomPoisRnd(n_sat*p_sim, 1);
            else
                N_fail = CustomPoisRnd(Non_k*p_sim, 1);
            end
        end
        
        % Update the number of available stock
        Non_k = max([Non_k - N_fail, 0]);
        
        %%% 2. Check Direct Resupply Arrival
        if cnt_lv == 0 % Arrive at this step
            % Update Non and Xq
            Non_k = Non_k + Q2;
            Xq2(Non_k+1,cnt_p,iter) = Xq2(Non_k+1,cnt_p,iter) + 1; % +1 for index
            
            % Update LV Parameters
            cnt_lv = -1;
            lv_cnt = lv_cnt + 1;
        elseif cnt_lv > 0 % Wait for arrival
            cnt_lv = cnt_lv - 1;
        end
        
        % 3. Check Q1 Resupply
        if cnt_p == 1
            % Demand
            n_Req = R1 + 1 - Non_k;
            if n_Req > 0
                n_dmd = ceil(n_Req/Q1);
                Xdmd(n_dmd+1,iter) = Xdmd(n_dmd+1,iter) + 1; % +1 for index
            else
                n_dmd = 0;
                Xdmd(1,iter) = Xdmd(1,iter) + 1;
            end
            
            % Update Xr
            Xr1(Non_k+1,iter) = Xr1(Non_k+1,iter) + 1; % +1 for index
            
            % Sample the number of available parking spares
            n_park = find(Pav_sum - rand >= 0, 1) - 1;
            Xav(n_park+1,iter) = Xav(n_park+1,iter) + 1; % +1 for index
            
            % Apply maximum feasible reorder #
            n_feas = min(n_park, n_dmd);
            Non_k = Non_k + n_feas*Q1;
            Xq1(Non_k+1,iter) = Xq1(Non_k+1,iter) + 1; % +1 for index
        end
        
        % 4. Check Q2 Resupply
        if cnt_lv == -1 && Non_k <= R2
            dT_LV = T0_LV + CustomExpRnd(mu_LV,1);
            dT_LV = ceil(dT_LV/dt_sim);
            if dTlv_max < dT_LV
                Xlv(end,iter) = Xlv(end,iter) + 1;
            else
                Xlv(dT_LV,iter) = Xlv(dT_LV,iter) + 1;
            end
            
            % Update Xr
            Xr2(Non_k+1,cnt_p,iter) = Xr2(Non_k+1,cnt_p,iter) + 1; % +1 for index
            
            % Save Remaining time step before arrival
            cnt_lv = dT_LV - 1;
        end
        
        % Save Stock Profile at current step after replinishment
        Non(k,iter) = Non_k;
        Xon(Non_k+1,iter) = Xon(Non_k+1,iter) + 1; % +1 for index
    end
end

%% Plot Simulation Result
xx_edge = -0.5:1:(Xmax+0.5);
dT_edge = 0.5:1:(dTlv_max+0.5);
dmd_edge = -0.5:1:(Dmax+0.5);

figure(1)
plot(time_sim, Non)

figure(2); hold on
histogram(Non(:),'Normalization','probability')
xlabel('Number of stock for entire period')
ylabel('Probability')

figure(3); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xq1,2), 'Normalization','probability')
xlabel('Number of stock right after Q1 resupply')
ylabel('Probability')

figure(4); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xr1,2), 'Normalization','probability')
xlabel('Number of stock at reordering of R1 resupply')
ylabel('Probability')

figure(5); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xq2,2), 'Normalization','probability')
xlabel('Number of stock right after Q2 resupply')
ylabel('Probability')

figure(6); hold on
histogram('BinEdges', xx_edge, 'BinCounts', sum(Xr2,2), 'Normalization','probability')
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


% mu_Q2 = idx_Q2/sum(idx_Q2);
% mu_R2 = idx_R2/sum(idx_R2);
% pi_Q2 = zeros(size(Ni_Q2));
% pi_R2 = zeros(size(Ni_R2));
% for i = 1:cnt_period
%     pi_Q2(:,i) = Ni_Q2(:,i)/sum(Ni_Q2(:,i));
%     pi_R2(:,i) = Ni_R2(:,i)/sum(Ni_R2(:,i));
% end

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
%.. State
Xmax = max(Q1+Q2+R1, Q1+Q2+R2);
x = Xmax:-1:0;
nx = length(x);

%.. Failure Transition Matrix
f = p_fail*dt_mc;
I = eye(nx);
Pf = zeros(nx);
for i = 1:nx
    if p_type == 0 % Constant Failure Rate
        Pf(i:end,i) = poisspdf(0:nx-i, n_sat*f)';        
    else % Stock Level Dependant Failure Rate
        if (nx-i) > n_sat
            Pf(i:end,i) = poisspdf(0:nx-i, n_sat*f)';        
        else
            Pf(i:end,i) = poisspdf(0:nx-i, (nx-i)*f)';  
        end
    end
    Pf(end,i) = 1 - sum(Pf(1:end-1, i));
end

%.. InDirect Resupply Transition Matrix
nx1 = Q1+R1+1;
mx1 = ceil(nx1/Q1);
Pq1_m = zeros(mx1*Q1);
for i = 1:mx1 % row
    for j = i:mx1 % column
        x_idx = (Q1*(i-1) + 1):(Q1*i);
        y_idx = (Q1*(j-1) + 1):(Q1*j);

        %.. Prob
        if i == 1
            % Only single kappa
            p_ij = Kappa(j);
        else
            % Kappa_i - Kappa_i+1
            i0 = j-i+1;
            p_ij = Kappa(i0) - Kappa(i0 + 1);
        end

        Pq1_m(x_idx,y_idx) = p_ij*eye(Q1);
    end
end
r = Xmax - Q1 - R1;
Pq1 = [eye(r), zeros(r,mx1*Q1); zeros(mx1*Q1,r), Pq1_m];
Pq1 = Pq1(1:nx,1:nx);

%.. Direct Resupply Transition Matrix 
Pq2 = zeros(nx);
Pq2(1:Q2, 1:Q2) = eye(Q2);
Pq2(1:nx-Q2, Q2+1:nx) = eye(nx-Q2);

%.. Selection Matrix
Cm = zeros(nx);
Cm(nx-R2:end,nx-R2:end) = eye(R2+1);
Cp = I - Cm;

%.. Full Transition
lam = 1/mu_LV;
rho = (1 - exp(-lam*dt_sim));

%.. Including Q1 moment
A11 = Cp*Pq1*Pf;
A12 = rho*Cp*Pq1*Pq2*Pf;
A21 = Cm*Pq1*Pf;
A22 = (1-rho)*Pq1*Pf + rho*Cm*Pq1*Pq2*Pf;
A0 = [A11, A12; A21, A22];

%.. Non-Q1 duration 
A11 = Cp*Pf;
A12 = rho*Cp*Pq2*Pf;
A21 = Cm*Pf;
A22 = (1-rho)*Pf + rho*Cm*Pq2*Pf;
A1 = [A11, A12; A21, A22];

%.. One Cycle
P_AA = A0*A1^(cnt_plane-1) ;
pi_AA = limitdist(P_AA');

%..
pi_f = zeros(nx*2, cnt_plane);
pi_f(:,1) = pi_AA;
for k = 2:cnt_plane-1
    pi_f(:,k) = A1 * pi_f(:,k-1);
end
pi_f(:,end) = A0*pi_f(:,end-1);

pi_tt = (pi_f(1:nx,end) + pi_f(nx+1:end,end));
for k = 1:cnt_plane-1
    pi_tt = pi_tt + pi_f(1:nx,k) + pi_f(nx+1:end,k);
end
pi_tt = pi_tt/cnt_plane;

figure(2)
plot(x, pi_tt, 'g*')
% xlim([x(16)-0.5, x(1)+0.5])
legend('Sim.', 'Sol.')

A11 = Pf*Cp + rho*Pf*Pq2*Cm;
A12 = rho*Pf*Pq2;
A21 = (1-rho)*Pf*Cm;
A22 = (1-rho)*Pf;
AA = [A11, A12; A21, A22];
P_AA = AA^cnt_plane * blkdiag(Pq1,Pq1);
pi_AA = limitdist(P_AA');

pi_f = zeros(nx*2, cnt_plane);
pi_f(:,1) = pi_AA;
pi_f(:,2) = AA * blkdiag(Pq1,Pq1) * pi_AA;
for k = 3:cnt_plane
    pi_f(:,k) = AA * pi_f(:,k-1);
end

pi_r = zeros(nx, cnt_plane);
pi_r(:,1) = Cm*Pq1*pi_f(1:nx,1);
for k = 2:cnt_plane
    pi_r(:,k) = Cm*pi_f(1:nx,k);
end
pi_R2A = sum(pi_r,2)/sum(sum(pi_r,2));

pi_tt = Pq1*(pi_f(1:nx,1) + pi_f(nx+1:end,1));
for k = 2:cnt_plane
    pi_tt = pi_tt + pi_f(1:nx,k) + pi_f(nx+1:end,k);
end
pi_tt = pi_tt/cnt_plane;

pi_q1_before = pi_f(1:nx,1) + pi_f(nx+1:end,1);
pi_q1_after = Pq1*pi_q1_before;

figure(2)
plot(x, pi_tt, 'r*')
% xlim([x(16)-0.5, x(1)+0.5])
legend('Sim.', 'Sol.')

% figure(3)
% plot(x, pi_q1_after, 'r*')
% xlim([x(16)-0.5, x(1)+0.5])
% legend('Sim.', 'Sol.')
% 
% figure(4)
% plot(x, pi_q1_before, 'r*')
% xlim([x(16)-0.5, x(1)+0.5])
% legend('Sim.', 'Sol.')

% figure(5)
% plot(x, sum(pi_Q2a,2), 'r*')
% 
% figure(6)
% plot(x, pi_R2A, 'r*')

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