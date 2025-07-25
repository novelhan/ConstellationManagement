% Def
% Compute Exact State Distrubution of Single Source Direct Resupply Method
%
% Input
% f: failure rate per unit time step
% f_type: failure type
% Q: Reorder quantity
% R: Reorder level
% mu: inverse of mean of lead time dist.
% dt_mc: unit time step
% dt_lv: process time of lead time
% n_sat: # of nominal satellite
%
% Output
% x: State vector of Markov Chain
% PI: Set of Stationary State Distribution
% T: Set of time duration for cycles
%
% Reference
% Analysis and Design of Satellite Constellation Spare Strategy Using Markov Chain
% https://doi.org/10.48550/arXiv.2408.09250

function [x, PI, T] = ExactDirectEventBased(f, f_type, Q, R, mu, dt_mc, dt_lv, n_sat)
    %%% Step 1: Initialize the parameters
    %.. State
    xmax = Q + R; % Max State Level: bar(N_sat), p7
    x = (xmax:-1:0)';
    nx = length(x); % Dimension of state distribution: bar(N_sat)+1
    
    %.. Constant lead time
    m = floor(dt_lv/dt_mc); % Below Eq.21
    
    %%% Step 2: Compute pi_q and pi_r
    %.. Failure Transition Matrix, Eqs.11~12
    f = f*dt_mc; % Failure rate for T_mc period
    I = eye(nx); % Dim: (bar(N_sat)+1) x (bar(N_sat) + 1)
    Pf = zeros(nx); % Dim: (bar(N_sat)+1) x (bar(N_sat) + 1)
    for i = 1:nx
        if f_type == 0 % Constant Failure Rate
            Pf(i:end,i) = CustomPoisPdf(0:nx-i, n_sat*f)';

        else % Stock Level Dependant Failure Rate
            if (nx-i) > n_sat
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, n_sat*f)';        
            else
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, (nx-i)*f)';  
            end
        end
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end
    
    %.. Resupply Transition Matrix, Eq.18
    Pq = [ [eye(Q); zeros(R+1, Q)], [eye(R+1); zeros(Q,R+1)] ];

    %.. Selection Matrix, Eq.7
    Cp = zeros(nx);
    Cm = zeros(nx);
    Cp(1:Q,1:Q) = eye(Q);
    Cm(Q+1:end, Q+1:end) = eye(R+1);
    
    %.. Transition Matrix based on the Lead Time model
    %.. TMP Mat
    P1 = inv(I - Pf*Cp);
    P2 = inv(I - exp(-mu*dt_mc)*Pf);
    P3 = Pf^(m+1);

    %.. Lead Time Prob, Eq.15
    rho0 = 1 - exp(-mu*dt_mc);

    %.. Full Transition
    Prq = Cm*P1; % Eq.17 
    Pqr = rho0*Pq*P3*P2; % Eq.21
    Pqq = Pqr*Prq; % Eq.22

    %.. Conditional Dist.
    pi_q = limitdist(Pqq); % Prob. Dist after reorder arrives
    pi_r = Prq*pi_q; % Prob. Dist when reorder is made

    %%% Step 3: Compute pi_np and pi_wp
    %.. Weighted Dist. 
    pi_np = Cp*P1*pi_q; % Eq.23
    A1 = I;
    for i = 1:m
        A1 = A1 + Pf^i;
    end
    A2 = (1 - rho0)*P3*P2;
    pi_wp = (A1 + A2)*pi_r; % Eq.27, If m = 0, A1 + A2 = P2
    
    %.. Avg. Time
    T_np = sum(pi_np); % Eq.24 in unit of [T_mc]
    T_wp = sum(pi_wp); % Eq.26
    T_dr = T_np + T_wp; % Eq.29
    
    %.. Avg. Dist.
    pi_dr = (pi_np + pi_wp)/T_dr; % Eq.28
    pi_np = pi_np/T_np;
    pi_wp = pi_wp/T_wp;
    
    %.. Output
    PI.pi_dr = pi_dr;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    
    T.T_dr = T_dr*dt_mc; % Eq.24 in unit of [day]
    T.T_wp = T_wp*dt_mc; % Eq.26
    T.T_np = T_np*dt_mc; % Eq.29
end