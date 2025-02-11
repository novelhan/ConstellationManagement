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
    k_p = round(dt_park/dt_mc); % for each different period
    k_lv = round(dt_bias_lv/dt_mc);

    %.. Constant lead time and residual
    mu = 1/dt_mu_lv; % mu = 1/(mean)
    m_p = floor(dt_bias_lv/dt_park); % Eq.43
    c_a = k_lv - m_p*k_p; % c_a = T_lp/T_mc, Eq.43
    c_b = (m_p+1)*k_p - k_lv; % c_b = T_rp/T_mc, Eq.43
    
    %%% Step 2: Compute pi_q and pi_r
    %.. Failure Transition Matrix of Parking Orbit, Eq.41
    I = eye(nx);
    Pf = zeros(nx,nx);
    for i = 1:nx
        Pf(i:end,i) = Eta(1:nx-i+1)';
        Pf(end,i) = 1 - sum(Pf(1:end-1,i));
    end
    
    %.. Resupply Transition Matrix
    Pq = zeros(nx,nx);
    for k = 1:n_plane
        Pq(:,:) = [ [eye(Q); zeros(R+1, Q)], [eye(R+1); zeros(Q,R+1)] ];
    end
    
    %.. Selection Matrix
    Cm = blockdiag(zeros(nx-R-1),eye(R+1));
    Cp = eye(nx) - Cm;
    
    %.. Stepwise full transition matrix
    nd = 2 + (k_lv - 1); % # of state vector consisting the full state
    nxx = nd*nx; % Full state dimension
    Pd = zeros(nxx,nxx); % Full state transition matrix
    
    for i = 1:nd
        if i == 1
            Pd(1:nx,1:nx) = Cp;
            Pd(1:nx,nx+1:2*nx) = rho*Cp*Pq;
        elseif i == 2
            Pd(nx+1:2*nx,nx+1:2*nx,k) = (1-rho)*I;
            Pd(nx+1:2*nx,end-(nx-1):end) = I;
        elseif i == 3
            Pd(2*nx+1:3*nx,1:nx,k) = Cm;
        else
            Pd((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = I;
        end
    end
    
    %.. Segmentwise full transition matrix
    Pseg = Pd^k_p;
    
    %.. Minimum period transition matrix
    Ptr = Pseg*kron(eye(nd),Pf);
    
    %.. Stationary solution
    PI_seg = limitdist(Ptr);
    PI_full = zeros(nxx,k_p);
    pi_full = zeros(nx,k_p);
    PI_full(:,1) = kron(eye(nd),Pf)*PI_seg;
    pi_full(:,1) = sum(reshape(PI_full(:,1),nx,[]),2);
    for i = 2:k_p
        PI_full(:,i) = Pd*PI_full(:,i-1);
        pi_full(:,i) = sum(reshape(PI_full(:,i),nx,[]),2);
    end

    %.. Distribution of waiting period
    pi_wp = zeros(nx,k_p);
    for i = 1:k_p
        pi_wp(:,i) = sum(reshape(PI_full(nx+1:2*nx,i),nx,[]),2);
        pi_wp(:,i) = pi_wp(:,i)/sum(pi_wp(:,i));
    end

    %.. Distribution of non-waiting period
    pi_np = zeros(nx,k_p);
    for i = 1:k_p
        pi_np(:,i) = sum(reshape(PI_full(1:nx,i),nx,[]),2);
        pi_np(:,i) = pi_np(:,i)/sum(pi_np(:,i));
    end
    
    %.. Avg. Dist.
    pi_ir = mean(pi_full,2);

    %.. Parking Available Distribution, Eq.34
    Pdav = zeros(1,nx);
    for i = 1:nx
        % Probability of having stock level larger than (i-1)
        Pdav(i) = sum(pi_ir(1:nx+1-i));
    end

    % Prob. Dist after reorder arrives
    pi_q = zeros(nx,k_p);
    pi_q(:,1) = Pq*pi_wp(:,end);
    pi_q(:,1) = pi_q(:,1)/sum(pi_q(:,1));
    for i = 2:k_p
        
        pi_q(:,i) = Pq*pi_wp(:,i-1);
        pi_q(:,i) = pi_q(:,i)/sum(pi_q(:,i));
    end

    % Prob. Dist when reorder is made
    pi_r = zeros(nx,k_p);
    pi_r(:,1) = Cm*pi_np(:,end);
    pi_r(:,1) = pi_r(:,1)/sum(pi_r(:,1));
    for i = 2:k_p
        pi_r(:,i) = Cm*pi_np(:,i-1);
        pi_r(:,i) = pi_r(:,i)/sum(pi_r(:,i));
    end
    
    % Should be same as
    for i = 1:k_p
        pi_r(:,i) = PI_full(2*nx+1:3*nx,i);
        pi_r(:,i) = pi_r(:,i)/sum(pi_r(:,i));
    end

    %.. Output
    PI.pi_ir = pi_ir;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.Pdav = Pdav;
    T = 0;
end

%%
function p = limitdist(P,method)

%.. Select computation method
if nargin == 1
    method = 1;
end
    
if method == 2
    % State Space Reduction Method
    % Obtain the stationary probability distribution vector p of an irreducible, recurrent Markov
    % chain by state reduction. P is the transition probabilities matrix of a discrete-time Markov
    % chain or the generator matrix Q.
    % https://www.math.wustl.edu/~feres/Math450Lect04.pdf
    
    P = P'; % This code assumes the row vector notation.
    [ns, ~]=size(P);
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
    
else
    % Simple solution approach for linear equation
    % Rmk: Based on the several test, this approach has better computational speed w/ same accuracy
    [n, ~] = size(P);
    A = [eye(n)-P; ones(1,n)];
    b = [zeros(n,1); 1];
    p = A\b;
end

end