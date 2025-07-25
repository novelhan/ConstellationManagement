% Def
% Compute Exact State Distrubution of Parking Orbit Spares under InDirect Resupply Method for symmetric orbits
%
% Input
% Eta: Demand from each different In-Plane orbit at every contact period
% Q: Reorder quantity
% R: Reorder level
% dt_mc: unit time step
% dt_park: review period
% dt_mu_lv: exponential mean time of lead time
% dt_bias_lv: process time of lead time
%
% Output
% x: State vector of Markov Chain
% PI: Set of Stationary State Distribution
% T: Set of time duration for cycles
%
% Reference

function [x, PI, T] = ExactInDirectParkRatio(Eta, Q, R, dt_mc, dt_park, dt_mu_lv, dt_bias_lv)
    %%% Step 1: Initialize the parameters
    %.. State
    xmax = Q + R; % Max State Level: bar(N_sat), p.10
    x = (xmax:-1:0)';
    nx = length(x); % Dimension of state distribution: bar(N_sat)+1
    
    %.. Step Count
    k_p = round(dt_park/dt_mc);
    k_lv = round(dt_bias_lv/dt_mc);
    rho = (1 - exp(-dt_mc/dt_mu_lv));

    %%% Step 2: Compute pi_q and pi_r
    %.. Failure Transition Matrix of Parking Orbit, Eq.41
    eta = Eta;
    if length(eta) <= nx
        Eta = zeros(nx,1);
        Eta(1:length(eta)) = eta;
    else
%         Eta = eta(1:nx);
%         Eta(end) = Eta(end) + sum(eta(nx+1:end));
    end
    Eta = Eta/sum(Eta);
    
    I = eye(nx);
    Pf = zeros(nx,nx);
    for i = 1:nx
        Pf(i:end,i) = Eta(1:nx-i+1)';
        Pf(end,i) = 1 - sum(Pf(1:end-1,i));
    end
    
    %.. Resupply Transition Matrix
    Pq = [ [eye(Q); zeros(R+1, Q)], [eye(R+1); zeros(Q,R+1)] ];
    
    %.. Selection Matrix
    Cm = blkdiag(zeros(nx-R-1),eye(R+1));
    Cp = I - Cm;
    
    %.. Stepwise full transition matrix
    %.. Assume resupply -> distribution -> check reorder (for different simulation, order may be needed to be changed)
    nd = 2 + (k_lv - 1); % # of state vector consisting the full state
    nxx = nd*nx; % Full state dimension
    Pd = zeros(nxx,nxx); % Full state transition matrix
    Pd0 = zeros(nxx,nxx); % State Transition at the contact moment
    
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
    a_np = sum(mean(pi_np,2));
    a_wp = sum(mean(pi_wp,2));
    a_tp = sum(mean(pi_tp,2))/(k_lv-1)*k_lv; % Only k_lv - 1 terms are considered for the states

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
    pi_np = sum(pi_np,2)/sum(sum(pi_np,2));
    pi_wp = sum(pi_wp+pi_tp,2)/sum(sum(pi_wp+pi_tp,2));
    
    %.. Output
    PI.pi_full = pi_full;
    PI.pi_ir = pi_ir;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.Pdav = Pdav;
    
    %.. Compuite Avg. period of each conditional distribution.
    T_tp = k_lv*dt_mc;
    T_np = a_np/a_tp*T_tp;
    T_wp = a_wp/a_tp*T_tp;
    T_ir = T_np + T_wp + T_tp;
    
    T.T_ir = T_ir; 
    T.T_wp = T_wp + T_tp; 
    T.T_np = T_np; 
end
