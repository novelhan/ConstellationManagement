% Def
% Compute Exact State Distrubution of Parking Orbit Spares under InDirect Resupply Method
%
% Input
% Eta: Demand from In-Plane orbit at every contact period
% Q: Reorder quantity
% R: Reorder level
% dt_mc: unit time step
% dt_park: review period (duration unitl meet subsequent in-plane orbit)
% dt_mu_lv: exponential mean time of lead time
% dt_bias_lv: process time of lead time
%
% Output
% x: State vector of Markov Chain
% PI: Set of Stationary State Distribution
% T: Set of time duration for cycles
%
% Reference
% Analysis and Design of Satellite Constellation Spare Strategy Using Markov Chain
% https://doi.org/10.48550/arXiv.2408.09250

function [x, PI, T] = ExactInDirectParkTime(Eta, Q, R, dt_mc, dt_park, dt_mu_lv, dt_bias_lv)
    %%% Step 1: Initialize the parameters
    %.. State
    xmax = Q + R; % Max State Level: bar(N_sat), p.10
    x = (xmax:-1:0)';
    nx = length(x); % Dimension of state distribution: bar(N_sat)+1
    
    %.. Step Count
    k_p = round(dt_park/dt_mc); % Eq.40
    k_lv = round(dt_bias_lv/dt_mc);

    %.. Constant lead time and residual
    mu = 1/dt_mu_lv; % mu = 1/(mean)
    m_p = floor(dt_bias_lv/dt_park); % Eq.43
    c_a = k_lv - m_p*k_p; % c_a = T_lp/T_mc, Eq.43
    c_b = (m_p+1)*k_p - k_lv; % c_b = T_rp/T_mc, Eq.43
    
    %%% Step 2: Compute pi_q and pi_r
    %.. Failure Transition Matrix of Parking Orbit, Eq.41
    eta = Eta;
    if length(eta) <= nx
        Eta = zeros(nx,1);
        Eta(1:length(eta)) = eta;
    else
        Eta = eta(1:nx);
        Eta(end) = Eta(end) + sum(eta(nx+1:end));
    end
    Eta = Eta/sum(Eta);

    I = eye(nx);
    Pf = zeros(nx);
    for i = 1:nx
        Pf(i:end,i) = Eta(1:nx-i+1)';
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end

    %.. Resupply Transition Matrix, Eq.18
    Pq = [ [eye(Q); zeros(R+1, Q)], [eye(R+1); zeros(Q,R+1)] ];
    
    %.. Selection Matrix, Eq.7
    Cp = zeros(nx);
    Cm = zeros(nx);
    Cp(1:Q,1:Q) = eye(Q);
    Cm(Q+1:end, Q+1:end) = eye(R+1);

    %.. TMP Mat
    P1 = inv(I - exp(-mu*k_p*dt_mc)*Pf);
    P2 = Pf^m_p;

    %.. Lead Time Prob
    rho_p0 = 1 - exp(-mu*c_b*dt_mc); % Eq.44 
    rho_p1 = exp(-mu*c_b*dt_mc)*(1-exp(-mu*k_p*dt_mc)); % Eq.44

    %.. Full Transition
    Prq = Cm*Pf*inv(I-Cp*Pf); % Eq.42
    Pqr = Pq*P2*(rho_p0*I + rho_p1*Pf*P1); % Eq.45
    Pqq = Pqr*Prq;
        
    %.. Conditional Dist.
    pi_q = limitdist(Pqq); % Prob. Dist after reorder arrives
    pi_r = Prq*pi_q; % Prob. Dist when reorder is made

    %%% Step 3: Compute pi_np and pi_wp
    %.. Full duration term of the first T_lv duration of waiting period
    A1 = zeros(nx);
    for i = 0:m_p-1
        A1 = A1 + Pf^i; % I + Pf + Pf^2 + ... + Pf^(m_p-1)
    end
    A1 = k_p*A1; % Eq.47-1

    %.. Partial duration term of the first T_lv duration of waiting period
    A2 = c_a*Pf^m_p; % Eq.47-2

    %.. Remaing partial duration and full duration terms of waiting period
    ec = exp(-mu*dt_mc);
    rho0_c = ec*(ec^c_b - 1)/(ec - 1);
    A3 = rho0_c*Pf^m_p; % Eq.47-3
    
    rho1_c = ec^(c_b+1)*(ec^k_p - 1)/(ec - 1);
    A4 = rho1_c*P2*Pf*P1; % Eq.47-4

    %.. Distribution of waiting period, Eq.46
    pi_wp = (A1 + A2 + A3 + A4)*pi_r;
    T_wp = sum(pi_wp); 

    %.. Partial duration and full duration terms of non-waiting period
    rho_set = (1-exp(-mu*dt_mc))*exp(-mu*dt_mc*(0:k_p-1)); % Eq.15: rho0, rho1,...
    rho_set = rho_set/(1-exp(-mu*k_p*dt_mc)); % Eq.52: bar(rho0), bar(rho1),...
    c_set = c_b - (1:k_p); % Eq.51
    c_set = c_set + (c_set < 0)*k_p; % Eq.51
    
    %.. Distribution of non-waiting period
    pi_np = sum(rho_set.*c_set)*pi_q + k_p*Cp*Pf*inv(I-Cp*Pf)*pi_q; % Eq.49
    T_np = sum(pi_np);

    %.. Avg. Dist.
    T_ir = T_np + T_wp;
    pi_ir = (pi_np + pi_wp)/T_ir; % Eq.53
    pi_np = pi_np/T_np;
    pi_wp = pi_wp/T_wp;

    %.. Parking Available Distribution, Eq.34
    Pdav = zeros(nx,1);
    for i = 1:nx
        % Probability of having stock level larger than (i-1)
        Pdav(i) = sum(pi_ir(1:nx+1-i));
    end

    %.. Output
    PI.pi_ir = pi_ir;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.Prq = Prq;
    PI.Pqr = Pqr;
    PI.Pdav = Pdav;
    
    T.T_ir = T_ir*dt_mc;
    T.T_wp = T_wp*dt_mc;
    T.T_np = T_np*dt_mc;
end