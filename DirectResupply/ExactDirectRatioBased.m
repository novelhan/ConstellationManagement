% Def
% Compute Exact State Distrubution of Single Source Direct Resupply Method using the memoryless
% property of the exponential distribution. 
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

function [x, PI, T] = ExactDirectRatioBased(f, f_type, Q, R, mu, dt_mc, dt_lv, n_sat)
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
   
    %.. Lead Time Prob, Eq.15
    rho = 1 - exp(-mu*dt_mc);

    %.. Full Transition
    Ptr = [Cp*Pf, rho*Cp*Pq*Pf; 
           Pf^(m-1)*Cm*Pf, (1-rho)*Pf + rho*Pf^(m-1)*Cm*Pq*Pf];
    
    %.. Stationary sol
    PI_dr = limitdist(Ptr);
       
    %.. Dist. during the non-reordering period
    pi_np = PI_dr(1:nx);
    a_np = sum(pi_np);
    
    %.. Dist. during the reorder-waitting period
    pi_wp = PI_dr(nx+1:2*nx);
    a_wp = sum(pi_wp);
    
    %.. Dist. during the processing period
    pi_tp = zeros(nx,m-1);
    pi_tp(:,1) = Cm*Pf*pi_np + rho*Cm*Pq*Pf*pi_wp;
    for i = 2:(m-1)
        pi_tp(:,i) = Pf*pi_tp(:,i-1);
    end
    a_tp = sum(pi_tp(:,1));
    
    %.. Avg distribution
    pi_dr = PI_dr(1:nx) + PI_dr(nx+1:2*nx) + sum(pi_tp,2);
    
    %.. Normalization
    pi_np = pi_np/a_np;
    pi_wp = pi_wp/a_wp;
    pi_tp = pi_tp/a_tp;
    pi_dr = pi_dr/sum(pi_dr);
    
    %.. Prob. Dist when reorder is made
    pi_r = pi_tp(:,1);
    
    %.. Prob. Dist after reorder arrives
    pi_q = Pq*Pf*pi_wp; % It is assumed that simulation is done as failure -> resupply
    
    %.. (!!) Due to different definition for pi_wp, following equation is used (will be fixed)
    pi_wp = pi_wp*a_wp + sum(pi_tp,2)*a_tp;
    pi_wp = pi_wp/sum(pi_wp);
    
    %.. Compuite Avg. period of each conditional distribution.
    T_tp = m*dt_mc;
    T_np = a_np/(a_tp*m)*T_tp;
    T_wp = a_wp/(a_tp*m)*T_tp;
    T_dr = T_np + T_wp + T_tp;
    
    %.. Output
    PI.pi_dr = pi_dr;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    
    T.T_dr = T_dr; 
    T.T_wp = T_wp + T_tp; 
    T.T_np = T_np; 
end