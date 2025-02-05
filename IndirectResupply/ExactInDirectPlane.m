% Def
% Compute Exact State Distrubution of InPlane Orbit Spares under InDirect Resupply Method
%
% Input
% f_mc: failure rate per unit time step
% f_type: failure type
% n_sat: # of nominal satellite
% kappa: Parking orbit spare availability probability distribution
% Q: Reorder quantity
% R: Reorder level
% dt_mc: unit time step
% dt_plane: review period (duration unitl meet subsequent parking orbit)
%
% Output
% x: State vector of Markov Chain
% PI: Set of Stationary State Distribution
% T: Set of time duration for cycles
%
% Reference
% Analysis and Design of Satellite Constellation Spare Strategy Using Markov Chain
% https://doi.org/10.48550/arXiv.2408.09250


function [x, PI, T] = ExactInDirectPlane(f_sim, f_type, n_sat, kappa, Q, R, dt_mc, dt_plane)
    %%% Step 1: Initialize the parameters
    %.. State
    xmax = Q + R; % Max State Level: bar(N_sat), p.10
    x = (xmax:-1:0)';
    nx = length(x); % Dimension of state distribution: bar(N_sat)+1
    

    %.. Constant reorder period
    k = round(dt_plane/dt_mc); % Eq.32

    %%% Step 2: Compute pi_q and pi_r
    %.. Failure Transition Matrix (Prq), Eqs.11~12
    f_mc = f_sim*dt_mc; % Failure rate per unit time step
    I = eye(nx);
    Pf = zeros(nx);
    for i = 1:nx
        if f_type == 0 % Constant Failure Rate
            Pf(i:end,i) = CustomPoisPdf(0:nx-i, n_sat*f_mc)';        
        else % Stock Level Dependant Failure Rate
            if (nx-i) > n_sat
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, n_sat*f_mc)';        
            else
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, (nx-i)*f_mc)';  
            end
        end
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end
    Prq = Pf^k; % Eq.33

    %.. Resupply Transition Matrix (Pqr), Eq. 36
    m = ceil(nx/Q); % This variable is needed only for implementation of the code
    Pq = zeros(m*Q);
    for i = 1:m % row
        for j = i:m % column
            x_idx = (Q*(i-1) + 1):(Q*i);
            y_idx = (Q*(j-1) + 1):(Q*j);
            
            %.. Prob
            if i == 1
                % Only single kappa
                p_ij = kappa(j);
            else
                % Kappa_i - Kappa_i+1
                i0 = j-i+1;
                p_ij = kappa(i0) - kappa(i0 + 1);
            end
            
            Pq(x_idx,y_idx) = p_ij*eye(Q);
        end
    end
    Pqr = Pq(1:nx,1:nx);

    %.. Full Transition
    Pqq = Pqr*Prq;
        
    %.. Conditional Dist. Eq.37
    pi_q = limitdist(Pqq); % Prob. Dist right before RAAN Contact
    pi_r = Prq*pi_q; % Prob. Dist right after RAAN Contact
    
    %%% Step 3: Compute pi_ir
    %.. Demand Distributio, Eq.39
    pi_bf = zeros(m*Q,1);
    pi_bf(1:nx) = pi_r; % Dim. m*Q x 1
    pi_dmd = (sum(reshape(pi_bf,Q,m)))'; % Dim. m x 1

    %.. Weighted Dist. Eq.38
    A1 = I;
    for i = 1:k-1
        A1 = A1 + Pf^i;
    end
    pi_ir = A1*pi_q/k;

    %.. Output
    PI.pi_ir = pi_ir;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.pi_dmd = pi_dmd;
    PI.Prq = Prq;
    PI.Pqr = Pqr;
    
    T.T_ir = k*dt_mc; % T_ir = T_plane = k*dt_mc
end