% Def
% Compute Exact State Distrubution of InPlane Orbit Spares under InDirect Resupply Method with
% asuumption of asymmetic parking orbits distribution.
%
% Input
% f_mc: failure rate per unit time step
% f_type: failure type
% n_sat: # of nominal satellite
% n_park: # of parking orbits
% kappa: Spare availability probability distribution of each different parking orbit
% Q: Reorder quantity
% R: Reorder level for each different parking orbits
% dt_mc: unit time step
% dt_plane: review period for each different parking orbits
%
% Output
% x: State vector of Markov Chain
% PI: Set of Stationary State Distribution
% T: Set of time duration for cycles
%
% Reference

function [x, PI, T] = AsymInDirectPlane(f_mc, f_type, n_sat, n_park, Kappa, Q, R, dt_mc, dt_plane)
    %.. Check basic dimension of the input
    if n_park ~= length(R)
        error('Dim. error on R')
    elseif n_park ~= length(dt_plane)
        error('Dim. error on dt_plane')
    end

    %.. State
    xmax = Q + max(R);
    x = (xmax:-1:0)';
    nx = length(x);

    %.. Constant reorder period
    c_plane = round(dt_plane/dt_mc);

    %.. Stepwise Failure Transition Matrix
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
    
    %.. Set of power of failure matrix
    c_max = max(c_plane);
    Pf_set = zeros(nx,nx,c_max); % Save each power for later usage
    Pf_set(:,:,1) = Pf;
    for i = 2:c_max
        Pf_set(:,:,i) = Pf*Pf_set(:,:,i-1);
    end
    Pf_sum = cumsum(Pf_set,3);
    
    %.. Failure Transition Matrix between the contact points
    %.. ex. Pq2r(:,:,1) = Transition from PIq(1) to PIr(2) (After #1 to Before #2)
    %.. ex. Pq2r(:,:,end) = Transition from PIq(n) to PIr(1) (After #n to Before #1)
    Pq2r = zeros(nx,nx,n_park);
    for i = 1:n_park
        Pq2r(:,:,i) = Pf_set(:,:,c_plane(i));
    end

    %.. Resupply Transition Matrix by each parking orbit (Pqr)
    %.. ex. Pr2q(:,:,1) = Replenishment matrix by #1 parking orbit (Before #1 to After #1)
    Pr2q = zeros(nx,nx,n_park);
    for k = 1:n_park % parking orbit
        nx_k = Q + R(k) + 1;
        m_k = ceil(nx_k/Q);
        Pq = zeros(m_k*Q,m_k*Q);
        for i = 1:m_k % row
        for j = i:m_k % column
            x_idx = (Q*(i-1) + 1):(Q*i);
            y_idx = (Q*(j-1) + 1):(Q*j);

           %.. Prob
            if i == 1
                % Only single kappa
                p_ij = Kappa(k,j);
            else
                % Kappa_i - Kappa_i+1
                i0 = j-i+1;
                p_ij = Kappa(k,i0) - Kappa(k,i0+1);
            end

            Pq(x_idx,y_idx) = p_ij*eye(Q); 
        end
        end
        
        % Match the dimension.
        Pq = Pq(1:nx_k,1:nx_k);
        if nx_k < nx
            Pq = blkdiag(eye(nx-nx_k), Pq); 
        end
        Pr2q(:,:,k) = Pq;
    end

    %.. Full Transition (Start from Before #1)
    Pcycle = I;
    for i = 1:n_park
        Pcycle = Pq2r(:,:,i)*Pr2q(:,:,i)*Pcycle;
    end
    
    %.. Conditional Dist.
    pi_r = zeros(nx,n_park); % Prob. Dist right before RAAN Contact
    pi_q = zeros(nx,n_park); % Prob. Dist right after RAAN Contact
    
    pi_r(:,1) = limitdist(Pcycle);
    pi_q(:,1) = Pr2q(:,:,1)*pi_r(:,1);
    for i = 2:n_park
        pi_r(:,i) = Pq2r(:,:,i-1)*pi_q(:,i-1); % After i-1 to Before i
        pi_q(:,i) = Pr2q(:,:,i)*pi_r(:,i); % Before i to After i
    end
    
    %.. Avg. Dist for Repeated Cycle (after each replenishment)
    pi_ir = zeros(nx,1);
    for i = 1:n_park
        pi_ir = pi_ir + (I + Pf_sum(:,:,c_plane(i)-1))*pi_q(:,i);
    end
    pi_ir = pi_ir/sum(c_plane); % Normalization;
    
    %.. Avg. Demand for each parking orbits
    pi_dmd = zeros(ceil(nx/Q),n_park);
    for i = 1:n_park
        nx_k = Q + R(i) + 1;
        m_k = ceil(nx_k/Q);
        pi_rb = zeros(m_k*Q,1); % elements of pi_r are grouped in batch Q
        pi_rb(1:nx_k) = pi_r(end-nx_k+1:end,i); % Dim. m_k*Q x 1
        pi_dmd(1:m_k,i) = (sum(reshape(pi_rb,Q,m_k)))'; % Dim. m x 1
        pi_dmd(1,i) =  pi_dmd(1,i) + sum(pi_r(1:(nx-nx_k),1));
%         pi_dmd(1,i) =  pi_dmd(1,i) + (1 - sum(pi_dmd(:,i))); % same as above
    end
    
    %.. Avg. period for a repated cycle
    T_ir = sum(c_plane)*dt_mc;
    
    %.. Output
    PI.pi_ir = pi_ir;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.pi_dmd = pi_dmd;
    T.T_ir = T_ir;
    
end