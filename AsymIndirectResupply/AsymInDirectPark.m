% Def
% Compute Exact State Distrubution of Parking Orbit Spares under InDirect Resupply Method with
% asuumption of asymmetic parking orbits distribution.
%
% Input
% Eta: Demand from each different In-Plane orbit at every contact period
% Q: Reorder quantity
% R: Reorder level for each different in-plane segment
% n_plane: # of asy in=plane orbits
% dt_mc: unit time step
% dt_park: review period for each different in-plane segment
% dt_mu_lv: exponential mean time of lead time
% dt_bias_lv: process time of lead time
%
% Output
% x: State vector of Markov Chain
% PI: Set of Stationary State Distribution
% T: Set of time duration for cycles
%
% Reference

function [x, PI, T] = AsymInDirectPark(Eta, Q, R, n_plane, dt_mc, dt_park, dt_mu_lv, dt_bias_lv)
    %.. Check basic dimension of the input
    if n_plane ~= length(R)
        error('Dim. error on R')
    elseif n_plane ~= length(dt_park)
        error('Dim. error on dt_plane')
    end

    %%% Step 1: Initialize the parameters
    %.. State
    xmax = Q + max(R);
    x = (xmax:-1:0)';
    nx = length(x); % Dimension of state distribution: bar(N_sat)+1
    
    %.. Step Count
    k_p = round(dt_park/dt_mc); % for each different period
    k_lv = round(dt_bias_lv/dt_mc);
    rho = (1 - exp(-dt_mc/dt_mu_lv));
    
    %%% Step 2: Compute pi_q and pi_r
    %.. Failure Transition Matrix of Parking Orbit, Eq.41
    I = eye(nx);
    Pf = zeros(nx,nx,n_plane); % for each in-plane orbit
    for k = 1:n_plane
    for i = 1:nx
        Pf(i:end,i,k) = Eta(k,1:nx-i+1)';
        Pf(end,i,k) = 1 - sum(Pf(1:end-1,i,k));
    end
    end
    
    %.. Selection Matrix (for each in-plane orbit)
    Cp = zeros(nx,nx,n_plane);
    Cm = zeros(nx,nx,n_plane);
    
    for k = 1:n_plane
        Cm(:,:,k) = blkdiag(zeros(nx-R(k)-1),eye(R(k)+1));
        Cp(:,:,k) = I - Cm(:,:,k);
    end
    
    %.. Resupply Transition Matrix (Only Pq for max(R) is needed)
    Pq = [ [eye(Q); zeros(max(R)+1, Q)], [eye(max(R)+1); zeros(Q,max(R)+1)] ];
    
    %.. Stepwise full transition matrix (for each in-plane)
    % TODO: Rename variables and make a code readable
    nd = 2 + (k_lv - 1); % # of state vector consisting the full state
    nxx = nd*nx; % Full state dimension
    Pd = zeros(nxx,nxx,n_plane);
    Pd1 = zeros(nxx,nxx,n_plane); % Full state transition matrix for each in-plane
    
    for i = 1:nd
    for k = 1:n_plane
        if i == 1
            Pd1(1:nx,1:nx,k) = Cp(:,:,k);
            Pd1(1:nx,nx+1:2*nx,k) = rho*Cp(:,:,k)*Pq;
            
            Pd(1:nx,1:nx,k) = I;
            Pd(1:nx,nx+1:2*nx,k) = rho*Pq;
        elseif i == 2
            Pd1(nx+1:2*nx,nx+1:2*nx,k) = (1-rho)*I;
            Pd1(nx+1:2*nx,end-(nx-1):end,k) = I;
            
            Pd(nx+1:2*nx,nx+1:2*nx,k) = (1-rho)*I;
            Pd(nx+1:2*nx,end-(nx-1):end,k) = I;
            
        elseif i == 3
            Pd1(2*nx+1:3*nx,1:nx,k) = Cm(:,:,k);
            Pd1(2*nx+1:3*nx,nx+1:2*nx,k) = rho*Cm(:,:,k)*Pq;
        else
            Pd1((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx,k) = I;
            
            Pd((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx,k) = I;
        end
    end
    end
    
    %.. Segmentwise full transition matrix (for each in-plane)
    Pseg = zeros(nxx,nxx,n_plane);
    for k = 1:n_plane
        Pseg(:,:,k) = Pd(:,:,k)^(k_p(k)-1)*Pd1(:,:,k);
    end
    
    %.. Minimum period transition matrix
    Ptr = eye(nxx);
    
    for k = 1:n_plane
        Ptr = Pseg(:,:,k)*kron(eye(nd),Pf(:,:,k))*Ptr;
    end
    
    %.. Stationary solution
    k_full = sum(k_p);
    PI_full = zeros(nxx,k_full);
    pi_full = zeros(nx,k_full);
    
    cnt = 1;
    for k = 1:n_plane
        % Transition from k to k+1 part
        if k == 1
            PI_full(:,1) = Pd1(:,:,1)*kron(eye(nd),Pf(:,:,1))*limitdist(Ptr);
        else
            PI_full(:,cnt) = Pd1(:,:,k)*kron(eye(nd),Pf(:,:,k))*PI_full(:,cnt-1);
        end
        pi_full(:,cnt) = sum(reshape(PI_full(:,cnt),nx,[]),2);
        cnt = cnt + 1;
        
        for i = 2:k_p(k)
            PI_full(:,cnt) = Pd(:,:,k)*PI_full(:,cnt-1);
            pi_full(:,cnt) = sum(reshape(PI_full(:,cnt),nx,[]),2);
            cnt = cnt + 1;
        end 
    end
    
    %.. Distribution of waiting period
    pi_wp = zeros(nx,k_full);
    for i = 1:k_full
        pi_wp(:,i) = sum(reshape(PI_full(nx+1:2*nx,i),nx,[]),2);
        % pi_wp(:,i) = pi_wp(:,i)/sum(pi_wp(:,i));
    end
    
    %.. Distribution of non-waiting period
    pi_np = zeros(nx,k_full);
    for i = 1:k_full
        pi_np(:,i) = sum(reshape(PI_full(1:nx,i),nx,[]),2);
        % pi_np(:,i) = pi_np(:,i)/sum(pi_np(:,i));
    end
        
    %.. Avg. Dist.
    pi_ir = mean(pi_full,2);

    %.. Parking Available Distribution, Eq.34
    Pdav = zeros(nx,n_plane);
    tmp = cumsum(k_p);
    for k = 1:n_plane
    for i = 1:nx
        % Probability of having stock level larger than (i-1)
        Pdav(i,k) = sum(pi_full(1:nx+1-i,tmp(k)));
    end
    end
    
    % Prob. Dist after reorder arrives
    pi_q = zeros(nx,k_full);
    cnt = 1;
    for k = 1:n_plane
        for i = 1:k_p(k)
            if cnt == 1
                pi_q(:,1) = Pq*pi_wp(:,end);
            else
                pi_q(:,cnt) = Pq*pi_wp(:,cnt-1);
            end
            % pi_q(:,cnt) = pi_q(:,cnt)/sum(pi_q(:,cnt));
            cnt = cnt + 1;
        end
    end

    % Prob. Dist when reorder is made
    pi_r = zeros(nx,n_plane);
    for k = 1:n_plane
        if k == 1
            pi_r(:,k) = PI_full(2*nx+1:3*nx,1);
        else
            pi_r(:,k) = PI_full(2*nx+1:3*nx,1+tmp(k-1));
        end
%         pi_r(:,k) = pi_r(:,k)/sum(pi_r(:,k));
    end

    %.. Output
    PI.PI_full = PI_full;
    PI.pi_full = pi_full;
    PI.pi_ir = pi_ir;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.Pdav = Pdav;
    
    T = 0;
end