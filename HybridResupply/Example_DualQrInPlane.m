close all
clear all
clc
% T 구간에서 T/2 를 2번가는 것과, T를 한번가는것 그리고 T/n을 n번 가는것 같은지 확인
% 포아송, 무한스테이트면 같을 것 같은데 한정스테이트, 스테이트 의존 고장모델도 같은지 확인 아마 아닐듯
%% Test Param
%.. Sim time
dt_sim      =   1;                  % [day]
time_sim    =   0:dt_sim:365*50000;

%.. Marcov Chain Period
dt_mc       =   dt_sim;             % [day]

%.. In-plane Period (Time duration of in-plane for subsequent parking contact, RAAN Drift time)
dt_plane   =   40;                 % [day]
cnt_plane  =   round(dt_plane/dt_sim);

%.. Failure rate
p_fail      =   0.55/365;           % [#/day]
p_sim       =   p_fail * dt_sim;    % [#/dt_sim]
p_mc        =   p_fail * dt_mc;     % [#/dt_mc]
p_type      =   1; % 0 for const, 1 for state dependant

%.. Hybrid In-plane (Q1,R1,Q2,R2) Policy Parameter
Q1  =   4;
R1  =   40;
Q2  =   2;
R2  =   38;
n_sat  =  40; % Nominal Satellite Count

%.. State Parameter
Xmax = Q1+Q2+R1; % Max State Level
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
dt_LV   =   20;
dTlv_max = ceil((dt_LV + mu_LV + 5*mu_LV)/dt_mc); % Max Lead Time Bin (mean + 5*sigma) [dt_sim]
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
            dT_LV = dt_LV + CustomExpRnd(mu_LV,1);
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

figure(7); hold on
histogram('BinEdges', dmd_edge, 'BinCounts', sum(Xdmd,2),'Normalization','probability')
xlabel('Number of demand')
ylabel('Probability')


Pbq2 = Xq2;
Pbr2 = Xr2;
for i = 1:cnt_plane
    Pbq2(:,i) = Pbq2(:,i)/sum(Pbq2(:,i));
    Pbr2(:,i) = Pbr2(:,i)/sum(Pbr2(:,i));
end

figure(8); hold on
for m = 1:(R2+1)
    plot(1:cnt_plane, Pbr2(m,:))
end

figure(9); hold on
for m = 1:(Xmax+1)
    plot(1:cnt_plane, [Pbq2(m,2:end), Pbq2(m,1)])
end

%% Run Analysis Code
ParaFail.dt_mc = dt_mc;
ParaFail.f_mc = p_mc;
ParaFail.f_type = p_type;
ParaFail.n_sat = n_sat;

%.. Direct Resupply Parameter
ParaDirect.mu = mu_LV;
ParaDirect.dt_lv = dt_LV;
ParaDirect.Q2 = Q2;
ParaDirect.R2 = R2;

%.. Indirect Resupply Parameter
ParaIndirect.dt_plane = dt_plane;
ParaIndirect.kappa = Kappa;
ParaIndirect.Q1 = Q1;
ParaIndirect.R1 = R1;

%..
ParaDim.n_cut = min(R1,R2);
% ParaDim.n_cut = 5;
% ParaDim.d_cut = dt_LV/dt_mc - 1;
ParaDim.d_cut = 3;
tic
[x, PI] = ExactDualInPlane(ParaFail, ParaDirect, ParaIndirect, ParaDim);
s = 'r*';

% [x, PI, Pf] = ExactDualInPlane1(ParaFail, ParaDirect, ParaIndirect, ParaDim);
% s = 'mx';
toc


figure(2); hold on
plot(x, sum(PI.pi_hr,2)/cnt_plane, s)
legend('Sim.', 'Sol.', 'location', 'best')

figure(3); hold on
plot(x, PI.pi_q1, s)
legend('Sim.', 'Sol.', 'location', 'best')

figure(4); hold on
plot(x, PI.pi_r1, s)
legend('Sim.', 'Sol.', 'location', 'best')

figure(5); hold on
plot(x, PI.pi_q2, s)
legend('Sim.', 'Sol.', 'location', 'best')

figure(6); hold on
plot(x, PI.pi_r2, s)
legend('Sim.', 'Sol.', 'location', 'best')

figure(7); hold on
plot(0:length(PI.pi_dmd)-1, PI.pi_dmd, s)
legend('Sim.', 'Sol.', 'location', 'best')

figure(8); hold on
set(gca,'ColorOrderIndex',1)
PI_R2 = flip(PI.PI_R2);
for m = 1:length(x)
    plot(1:cnt_plane, PI_R2(m,:), 'o')
end

PI_Q2 = flip(PI.PI_Q2);
figure(9); hold on
set(gca,'ColorOrderIndex',1)
for m = 1:length(x)
    plot(1:cnt_plane, PI_Q2(m,:), 'x')
end

% nx = length(x);
% figure(10); hold on
% set(gca,'ColorOrderIndex',1)
% for m = 1:length(x)
%     plot(1:cnt_plane, PI.PI_hr(m+nx,:), 'o')
% end
% 
% for m = 1:length(x)
%     plot(1:cnt_plane, PI.PI_hr(end-2*nx+m,:), 'x')
% end

% [ndim,ny] = size(PI.PI_hr);
% A = zeros(nx,ndim/nx,ny);
% for i = 1:ny
%     A(:,:,i) = reshape(PI.PI_hr(:,i),nx,[]);
%     for j = 1:ndim/nx
%         A(:,j,i) = A(:,j,i);
%     end
% end


