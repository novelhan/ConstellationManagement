%%
function [PI_i, PI_p, T_i, T_p, err] = ExactInDirectProb(iter_max, ParaInPlane, ParaParking)
    %.. Assume Perfect Availability for initialization
    nx_i = ParaInPlane.Q + ParaInPlane.R + 1;
    Pdav = ones(1,nx_i);
    PI_i = ExactInDirectPlane(Pdav, ParaInPlane);
    PI_p = ExactInDirectPark(PI_i.pi_dmd, ParaParking);

    %.. Fixed Point Root Finding
    xx_pre = [PI_i.pi_ir; PI_p.pi_ir];
    err = zeros(1,iter_max);
    for iter = 1:iter_max
        %.. Update the solution
        [PI_i, T_i] = ExactInDirectPlane(PI_p.Pdav, ParaInPlane);
        [PI_p, T_p] = ExactInDirectPark(PI_i.pi_dmd, ParaParking);
        
        %.. Check Convergence
        xx = [PI_i.pi_ir; PI_p.pi_ir];
        err(iter) = norm(xx - xx_pre);
        disp(['Iteration:', num2str(iter), ',  Error:', num2str(err(iter))])
        if err(iter) < 1e-5
            err = err(1:iter);
            break;
        else
            xx_pre = xx;
        end
    end
end

%% In-Plane Analysis
function [PI, T] = ExactInDirectPlane(Pavd, Para)
    %.. Parameter
    Q = Para.Q;
    R = Para.R;
    dt_mc = Para.dt_mc;
    dt_plane = Para.dt_plane;
    n_sat = Para.n_sat;

    %.. State
    xmax = Q + R;
    x = (xmax:-1:0)';
    nx = length(x);
    m = ceil(nx/Q);
    
    %.. Constant reorder period
    c_plane = round(dt_plane/dt_mc);
    
    %.. Failure Transition Matrix (Prq)
    f = Para.f_sim*dt_mc;
    I = eye(nx);
    Pf = zeros(nx);
    for i = 1:nx
        if Para.f_type == 0 
            % Constant Failure Rate
            Pf(i:end,i) = poisspdf(0:nx-i, n_sat*f)';        

        else                
            % Stock Level Dependant Failure Rate
            if (nx-i) > n_sat
                Pf(i:end,i) = poisspdf(0:nx-i, n_sat*f)';        
            else
                Pf(i:end,i) = poisspdf(0:nx-i, (nx-i)*f)';  
            end
        end
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end
    Prq = Pf^c_plane;

    %.. Resupply Transition Matrix (Pqr)
    Pq = zeros(m*Q);
    Pq1 = zeros(m*Q);
    for i = 1:m
        for j = i:m
            x_idx = (Q*(i-1) + 1):(Q*i);
            y_idx = (Q*(j-1) + 1):(Q*j);
            
            %.. Prob
            i0 = j-i+1;
            p_ij = Pavd(i0);
            for k = 1:i-1
                p_ij = p_ij*(1-Pavd(i0+k));
            end
            
            %.. Prob
            if i == 1
                % Only single kappa
                p_ij1 = Pavd(j);
            else
                % Kappa_i - Kappa_i+1
                i0 = j-i+1;
                p_ij1 = Pavd(i0) - Pavd(i0 + 1);
            end
            
            Pq(x_idx,y_idx) = p_ij*eye(Q);
            Pq1(x_idx,y_idx) = p_ij1*eye(Q);
        end
    end
    Pqr = Pq(1:nx,1:nx);
    Pqr1 = Pq1(1:nx,1:nx);

    %.. Full Transition
    Pqq = Pqr*Prq;
        
    %.. Conditional Dist.
    pi_q = limitdist(Pqq'); % Prob. Dist right before RAAN Contact
    pi_r = Prq*pi_q; % Prob. Dist right after RAAN Contact
    pi_bf = zeros(m*Q,1);
    pi_bf(1:nx) = pi_r;
    pi_dmd = (sum(reshape(pi_bf,Q,m)))';

    %.. Weighted Dist.
    Pir = I;
    for i = 1:c_plane-1
        Pir = Pir + Pf^i;
    end
    Pir = Pir/c_plane;
    pi_ir = Pir*pi_q;

    %.. Output
    PI.pi_ir = pi_ir;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.pi_dmd = pi_dmd;
    
    T.T_ir = c_plane*dt_mc;
end

%% Parking-Plane Analysis
function [PI, T] = ExactInDirectPark(Pdmd, Para)
    %%% Step 1: Initialize the parameters
    %.. Parameter
    Q = Para.Q;
    R = Para.R;
    dt_mc = Para.dt_mc;
    dt_park = Para.dt_park;
    dt_lv = Para.dt_lv;

    %.. State
    xmax = Q + R;
    x = (xmax:-1:0)';
    nx = length(x);
    
    %.. Step Count
    c_park = round(dt_park/dt_mc);
    c_lv = round(dt_lv/dt_mc);

    %.. Constant lead time and residual
    mu = 1/Para.mu_lv;
    m_p = floor(dt_lv/dt_park);
    c_a = c_lv - m_p*c_park;
    c_b = (m_p+1)*c_park - c_lv;

    %.. Failure Transition Matrix
    I = eye(nx);
    Pf = zeros(nx);
    eta = zeros(nx,1);
    eta(1:length(Pdmd)) = Pdmd;
    eta = eta/sum(eta);

    for i = 1:nx
        Pf(i:end,i) = eta(1:nx-i+1)';
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end

    %.. Resupply Transition Matrix
    Pq = [ [eye(Q); zeros(R+1, Q)], [eye(R+1); zeros(Q,R+1)] ];
    
    %%% Step 2: Compute pi_q and pi_r
    %.. Selection Matrix
    Cp = zeros(nx);
    Cm = zeros(nx);
    Cp(1:Q,1:Q) = eye(Q);
    Cm(Q+1:end, Q+1:end) = eye(R+1);

    %.. Lead Time Prob
    P1 = inv(I - exp(-mu*c_park*dt_mc)*Pf);
    P2 = Pf^m_p;

    rho_p0 = 1 - exp(-mu*c_b*dt_mc);  
    rho_p1 = exp(-mu*c_b*dt_mc)*(1-exp(-mu*c_park*dt_mc));

    %.. Full Transition
    Prq = Cm*Pf*inv(I-Cp*Pf);
    Pqr = Pq*P2*(rho_p0*I + rho_p1*Pf*P1);
    Pqq = Pqr*Prq;
        
    %.. Conditional Dist.
    pi_q = limitdist(Pqq'); % Prob. Dist after reorder arrives
    pi_r = Prq*pi_q; % Prob. Dist when reorder is made

    %%% Step 3: Compute pi_np and pi_wp
    %.. Full duration term of the first T_lv duration of waiting period
    A1 = zeros(nx);
    for i = 0:m_p-1
        A1 = A1 + Pf^i; % I + Pf + Pf^2 + ... + Pf^(m_p-1)
    end
    A1 = c_park*A1;

    %.. Partial duration term of the first T_lv duration of waiting period
    A2 = c_a*Pf^m_p;

    %.. Remaing partial duration and full duration terms of waiting period
    ec = exp(-mu*dt_mc);
    et = exp(-mu*c_park*dt_mc);
    rho0_c = ec*(ec^c_b - 1)/(ec - 1);
    rho1_c = ec^(c_b+1)*(ec^c_park - 1)/(ec - 1);
    A3 = rho0_c*Pf^m_p;
    A4 = rho1_c*P2*Pf*P1;

    %.. Distribution of waiting period
    pi_wp = (A1 + A2 + A3 + A4)*pi_r;
    T_wp = sum(pi_wp);

    %.. Partial duration and full duration terms of non-waiting period
    rho_set = (1-exp(-mu*dt_mc))*exp(-mu*dt_mc*(0:c_park-1)); %
    rho_set = rho_set/(1-exp(-mu*c_park*dt_mc));
    c_set = c_b - (1:c_park);
    c_set = c_set + (c_set < 0)*c_park;
    
    %.. Distribution of non-waiting period
    pi_np = sum(rho_set.*c_set)*pi_q + c_park*Cp*Pf*inv(I-Cp*Pf)*pi_q;
    T_np = sum(pi_np);

    %.. Avg. Dist.
    T_ir = T_np + T_wp;
    pi_ir = (pi_np + pi_wp)/T_ir;
    pi_np = pi_np/T_np;
    pi_wp = pi_wp/T_wp;

    %.. Parking Available Distribution
    Pdav = zeros(nx,1);
    for i = 1:nx
        Pdav(i) = sum(pi_ir(1:nx+1-i));
    end

    %.. Output
    PI.pi_ir = pi_ir;
    PI.pi_np = pi_np;
    PI.pi_wp = pi_wp;
    PI.pi_q = pi_q;
    PI.pi_r = pi_r;
    PI.Pdav = Pdav;
    
    T.T_ir = T_ir*dt_mc;
    T.T_wp = T_wp*dt_mc;
    T.T_np = T_np*dt_mc;
end

%%
function p = limitdist(P)
%Obtain the stationary probability distribution
%vector p of an irreducible, recurrent Markov
%chain by state reduction. P is the transition
%probabilities matrix of a discrete-time Markov
%chain or the generator matrix Q.
% https://www.math.wustl.edu/~feres/Math450Lect04.pdf

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
end