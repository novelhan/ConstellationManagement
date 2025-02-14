% Def
% Compute Exact State Distrubution of Single Source Direct Resupply Method
%
% Input
% ParaInPlane: In-Plane orbit parameter structure
%   f_ref: failure rate per unit time step
%   f_type: failure type
%   N_sat: # of nominal satellite
%   Q: Reorder quantity
%   R: Reorder level
%   dt_mc: unit time step
%   mu_lv: mean of lead time dist.
%   dt_lv: process time of lead time
%   method: 0: Event based, 1: Ratio based method
%
% Output
% PI: Set of Stationary State Distribution
%   X: State vector of Markov Chain
%   pi_dr: Prob. dist. for the entire cycle
%   pi_np: Prob. dist. for the non-reordering period
%   pi_wp: Prob. dist. for the waitting period
%   pi_q: Prob. dist. right after Q replenishment
%   pi_r: Prob. dist. when reorder is made (when Xi = r)
% T: Set of time duration for cycles
%   T_wp: Avg. duration for waitting period
%   T_np: Avg. duration for non-reordering period
%   T_dr: Avg. duration for full reorder cycle
%
% Reference
% Analysis and Design of Satellite Constellation Spare Strategy Using Markov Chain
% https://doi.org/10.48550/arXiv.2408.09250

function [PI, T] = SolveDirectProb(ParaInPlane)
    %%% Step 1: Initialize the parameters
    %.. State
    Q = ParaInPlane.Q;
    R = ParaInPlane.R;
    xmax = Q + R; % Max State Level: bar(N_sat), p7
    x = (xmax:-1:0)';
    nx = length(x); % Dimension of state distribution: bar(N_sat)+1
    
    %.. Constant lead time
    mu_lv = ParaInPlane.mu_lv;
    dt_lv = ParaInPlane.dt_lv;
    dt_mc = ParaInPlane.dt_mc;
    m = floor(dt_lv/dt_mc); % Below Eq.21
    
    %%% Step 2: Compute pi_q and pi_r
    f_ref = ParaInPlane.f_ref;
    f_type = ParaInPlane.f_type;
    N_sat = ParaInPlane.N_sat;
    
    %.. Failure Transition Matrix, Eqs.11~12
    f_mc = f_ref*dt_mc; % Failure rate for T_mc period
    I = eye(nx); % Dim: (bar(N_sat)+1) x (bar(N_sat) + 1)
    Pf = zeros(nx); % Dim: (bar(N_sat)+1) x (bar(N_sat) + 1)
    for i = 1:nx
        if f_type == 0 % Constant Failure Rate
            Pf(i:end,i) = CustomPoisPdf(0:nx-i, N_sat*f_mc)';

        else % Stock Level Dependant Failure Rate
            if (nx-i) > N_sat
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, N_sat*f_mc)';        
            else
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, (nx-i)*f_mc)';  
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
    
    %.. Solve the problem using two different methods
    if ParaInPlane.method == 0 % Event Based
        %.. Transition Matrix based on the Lead Time model
        %.. TMP Mat
        P1 = inv(I - Pf*Cp);
        P2 = inv(I - exp(-dt_mc/mu_lv)*Pf);
        P3 = Pf^(m+1);

        %.. Lead Time Prob, Eq.15
        rho0 = 1 - exp(-dt_mc/mu_lv);

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

        %.. Avg. Time expressed in unit step count
        T_np = sum(pi_np); % Eq.24 in unit of [T_mc]
        T_wp = sum(pi_wp); % Eq.26
        T_dr = T_np + T_wp; % Eq.29

        %.. Avg. Dist.
        pi_dr = (pi_np + pi_wp)/T_dr; % Eq.28
        pi_np = pi_np/T_np;
        pi_wp = pi_wp/T_wp;
        
        %.. Avg. Time expressed in day
        T_np = T_np*dt_mc;
        T_wp = T_wp*dt_mc;
        T_dr = T_dr*dt_mc;
        
    else % Ratio Based 
        %.. Lead Time Prob, Eq.15
        rho = 1 - exp(-dt_mc/mu_lv);

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
    end
    
    %.. Output
    PI.X = x;
    PI.pi_dr = pi_dr;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    
    T.T_dr = T_dr; % Eq.24 in unit of [day]
    T.T_wp = T_wp; % Eq.26
    T.T_np = T_np; % Eq.29
end