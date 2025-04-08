% Def
% Compute Exact State Distrubution of Parking Orbit Spares under InDirect Resupply Method
%
% Input
% Eta: Demand from In-Plane orbit at every contact period
% ParaParking: Parking orbit parameter structure
%   Q: Reorder quantity
%   R: Reorder level
%   dt_mc: unit time step
%   dt_park: review period (duration unitl meet subsequent in-plane orbit)
%   mu_lv: exponential mean time of lead time
%   dt_lv: process time of lead time
%   method: 0: Event based, 1: Ratio based method
%
% Output
% PI: Set of Stationary State Distribution
%   X: State vector of Markov Chain
%   pi_ir: Prob. dist. for the entire cycle
%   pi_np: Prob. dist. for the non-reordering period
%   pi_wp: Prob. dist. for the waitting period
%   pi_q: Prob. dist. right after Q replenishment
%   pi_r: Prob. dist. when reorder is made (when Xi = r)    
% T: Set of time duration for cycles
%   T_wp: Avg. duration for waitting period
%   T_np: Avg. duration for non-reordering period
%   T_ir: Avg. duration for full reorder cycle
%
% Reference
% Analysis and Design of Satellite Constellation Spare Strategy Using Markov Chain
% https://doi.org/10.48550/arXiv.2408.09250

function [PI, T] = SolveInDirectPark(Eta, ParaParking)
    %%% Step 1: Initialize the parameters
    %.. State
    Q = ParaParking.Q;
    R = ParaParking.R;
    xmax = Q + R; % Max State Level: bar(N_sat), p.10
    x = (xmax:-1:0)';
    nx = length(x); % Dimension of state distribution: bar(N_sat)+1
    
    %.. Step Count
    dt_mc = ParaParking.dt_mc;
    dt_park = ParaParking.dt_park;
    mu_lv = ParaParking.mu_lv;
    dt_lv = ParaParking.dt_lv;
    k_p = round(dt_park/dt_mc); % Eq.40
    k_lv = round(dt_lv/dt_mc);

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
    Eta = Eta;

    I = eye(nx);
    Pf = zeros(nx);
    for i = 1:nx
        Pf(i:end,i) = Eta(1:nx-i+1)';
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end

    %.. Resupply Transition Matrix, Eq.18
    Pq = [ [eye(Q); zeros(R+1, Q)], [eye(R+1); zeros(Q,R+1)] ];
    
    %.. Selection Matrix, Eq.7
    Cm = blkdiag(zeros(nx-R-1),eye(R+1));
    Cp = I - Cm;

    %.. Analysis Method Selection
    if isempty(ParaParking.method) || ParaParking.method == 0 % Event Based
%         if Q <= R
%             warning('Q > R for time based parking analysis')
%         end
            
        %.. Constant lead time and residual
        mu = 1/mu_lv; % mu = 1/(mean)
        m_p = floor(dt_lv/dt_park); % Eq.43
        c_a = k_lv - m_p*k_p; % c_a = T_lp/T_mc, Eq.43
        c_b = (m_p+1)*k_p - k_lv; % c_b = T_rp/T_mc, Eq.43
        
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
        T_avg = T_np + T_wp;
        pi_ir = (pi_np + pi_wp)/T_avg; % Eq.53
        pi_np = pi_np/T_np;
        pi_wp = pi_wp/T_wp;

        %.. Parking Available Distribution, Eq.34
        % Update papaer and clean the code
        pi_dmd_np = inv(I-Cp*Pf)*pi_q;
        pi_dmd_wp = (A1/k_p + (1-rho_p0)*Pf^m_p*P1)*pi_r;
        a = sum(pi_dmd_np);
        b = sum(pi_dmd_wp);
        pi_dmd = pi_dmd_np/(a+b) + pi_dmd_wp/(a+b);
        
        
        Pdav = zeros(nx,1);
        for i = 1:nx
            % Probability of having stock level larger than (i-1)
%             Pdav(i) = sum(pi_ir(1:nx+1-i));
            Pdav(i) = sum(pi_dmd(1:nx+1-i));
        end
        
        %.. Avg. Time expressed in day
        T_np = T_np*dt_mc;
        T_wp = T_wp*dt_mc;
        T_avg = T_avg*dt_mc;
        
    else % Ratio Based 
        %.. Stepwise full transition matrix
        %.. Assume resupply -> distribution -> check reorder (for different simulation, order may be needed to be changed)
        nd = 2 + (k_lv - 1); % # of state vector consisting the full state
        nxx = nd*nx; % Full state dimension
        Pd = zeros(nxx,nxx); % Full state transition matrix
        Pd0 = zeros(nxx,nxx); % State Transition at the contact moment
        rho = (1 - exp(-dt_mc/mu_lv)); % Prob. having replenishment within a single step.

        for i = 1:nd
            if i == 1
                Pd(1:nx,1:nx) = Cp;
                Pd(1:nx,nx+1:2*nx) = rho*Cp*Pq;

                Pd0(1:nx,1:nx) = Cp*Pf;
                Pd0(1:nx,nx+1:2*nx) = rho*Cp*Pf*Pq;
            elseif i == 2
                Pd(nx+1:2*nx,nx+1:2*nx) = (1-rho)*I;
                Pd(nx+1:2*nx,end-(nx-1):end) = I;

                Pd0(nx+1:2*nx,nx+1:2*nx) = (1-rho)*Pf;
                Pd0(nx+1:2*nx,end-(nx-1):end) = Pf;
            elseif i == 3
                Pd(2*nx+1:3*nx,1:nx) = Cm;
                Pd(2*nx+1:3*nx,nx+1:2*nx) = rho*Cm*Pq;

                Pd0(2*nx+1:3*nx,1:nx) = Cm*Pf;
                Pd0(2*nx+1:3*nx,nx+1:2*nx) = rho*Cm*Pf*Pq;
            else
                Pd((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = I;

                Pd0((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = Pf;
            end
        end
        %.. Minimum period transition matrix
        Ptr = Pd0*Pd^(k_p-1);
        
        %.. Stationary solution: idx=1 -> 1st step after review period
        PI_full = zeros(nxx,k_p);
        pi_full = zeros(nx,k_p);
        PI_full(:,1) = limitdist(Ptr);
        pi_full(:,1) = sum(reshape(PI_full(:,1),nx,[]),2);
        for i = 2:k_p
            PI_full(:,i) = Pd*PI_full(:,i-1);
            pi_full(:,i) = sum(reshape(PI_full(:,i),nx,[]),2);
        end

        %.. Distribution of non-waiting/waiting/processing period
        pi_np = zeros(nx,k_p);
        pi_wp = zeros(nx,k_p);
        pi_tp = zeros(nx,k_p);
        for i = 1:k_p
            pi_np(:,i) = PI_full(1:nx,i);
            pi_wp(:,i) = PI_full(nx+1:2*nx,i);
            pi_tp(:,i) = sum(reshape(PI_full(2*nx+1:end,i),nx,[]),2);
        end
        a_np = sum(sum(pi_np,2));
        a_wp = sum(sum(pi_wp,2));
        a_tp = sum(sum(pi_tp,2))/(k_lv-1)*k_lv; % Only k_lv - 1 terms are considered for the states
%         a_tp = sum(sum(pi_tp,2)); % Modifiy the state to use this code line

        %.. Avg. Dist.
        pi_ir = mean(pi_full,2);

        %.. Parking Available Distribution, Eq.34
        Pdav = zeros(nx,1);
        pi_dmd = pi_np(:,end-1) + (rho*Pq + (1-rho)*I )*pi_wp(:,end-1) + pi_tp(:,end-1);
        for i = 1:nx
            % Probability of having stock level larger than (i-1)
            Pdav(i) = sum(pi_dmd(1:nx+1-i,end));
        end

        % Prob. Dist after reorder arrives
        pi_q = zeros(nx,k_p);
        pi_q(:,1) = Pq*Pf*pi_wp(:,end);
        for i = 2:k_p
            pi_q(:,i) = Pq*pi_wp(:,i-1);
        end
        pi_q = sum(pi_q,2)/sum(sum(pi_q,2));

        % Prob. Dist when reorder is made
        pi_r = zeros(nx,k_p);
        for i = 1:k_p
            pi_r(:,i) = PI_full(2*nx+1:3*nx,i); 
        end
        pi_r = sum(pi_r,2)/sum(sum(pi_r,2));

        % Normalization
        pi_np = sum(pi_np,2)/a_np;
        pi_wp = sum(pi_wp+pi_tp,2)/(a_wp + a_tp);
        
        %.. Compuite Avg. period of each conditional distribution.
        T_tp = k_lv*dt_mc;
        T_np = a_np/a_tp*T_tp;
        T_wp = a_wp/a_tp*T_tp + T_tp;
        T_avg = T_np + T_wp;
    end
    
    %.. Output
    PI.X = x;
    PI.pi_avg = pi_ir;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.Pdav = Pdav;
    
    T.T_avg = T_avg;
    T.T_wp = T_wp;
    T.T_np = T_np;
end