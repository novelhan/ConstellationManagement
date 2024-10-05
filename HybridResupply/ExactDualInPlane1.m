function [x, PI, Pf] = ExactDualInPlane1(ParaFail, ParaDirect, ParaIndirect, ParaDim)
    %.. Failure Parameter
    dt_mc = ParaFail.dt_mc;
    f_mc = ParaFail.f_mc;
    f_type = ParaFail.f_type;
    n_sat = ParaFail.n_sat;
    
    %.. Direct Resupply Parameter
    mu = ParaDirect.mu;
    dt_lv = ParaDirect.dt_lv;
    Q2 = ParaDirect.Q2;
    R2 = ParaDirect.R2;
    
    %.. Indirect Resupply Parameter
    dt_plane = ParaIndirect.dt_plane;
    kappa = ParaIndirect.kappa;
    Q1 = ParaIndirect.Q1;
    R1 = ParaIndirect.R1;
    
    %.. Dimension Reduction Parameter
    n_cut = ParaDim.n_cut;
    d_cut = ParaDim.d_cut;
    
    %.. State
    xmax = Q1+Q2+R1;
    xmin = max(min(R1,R2) - n_cut,0);
    x = (xmax:-1:xmin)';
    nx = length(x);
    
    %.. Constant reorder period
    c_plane = round(dt_plane/dt_mc);
    Dmax = ceil(xmax/Q1);
    
    %.. LV Param
    lam = 1/mu;
    rho = (1 - exp(-lam*dt_mc));
    c_lv = round(dt_lv/dt_mc);
    
    %.. Failure Transition Matrix
    I = eye(nx);
    Pf = zeros(nx);
    for i = 1:nx
        if f_type == 0 % Constant Failure Rate
            Pf(i:end,i) = poisspdf(0:nx-i, n_sat*f_mc)';        
        else % Stock Level Dependant Failure Rate
            if (xmax+1-i) > n_sat
                Pf(i:end,i) = poisspdf(0:nx-i, n_sat*f_mc)';        
            else
                Pf(i:end,i) = poisspdf(0:nx-i, (xmax+1-i)*f_mc)';  
            end
        end
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end
    
    %.. InDirect Resupply Transition Matrix
    nx1 = Q1+R1+1;
    mx1 = ceil(nx1/Q1);
    Pq1_m = zeros(mx1*Q1);
    for i = 1:mx1 % row
        for j = i:mx1 % column
            x_idx = (Q1*(i-1) + 1):(Q1*i);
            y_idx = (Q1*(j-1) + 1):(Q1*j);

            %.. Prob
            if i == 1
                % Only single kappa
                p_ij = kappa(j);
            else
                % Kappa_i - Kappa_i+1
                i0 = j-i+1;
                p_ij = kappa(i0) - kappa(i0 + 1);
            end

            Pq1_m(x_idx,y_idx) = p_ij*eye(Q1);
        end
    end
    r = xmax - Q1 - R1;
    Pq1 = [eye(r), zeros(r,mx1*Q1); zeros(mx1*Q1,r), Pq1_m];
    Pq1 = Pq1(1:nx,1:nx);
    
    %.. Direct Resupply Transition Matrix 
    Pq2 = zeros(nx);
    Pq2(1:Q2, 1:Q2) = eye(Q2);
    Pq2(1:nx-Q2, Q2+1:nx) = eye(nx-Q2);
    
    %.. Selection Matrix
    Cp = zeros(nx);
    Cp(1:Q1+Q2+R1-R2,1:Q1+Q2+R1-R2) = eye(Q1+Q2+R1-R2);
    Cm = I - Cp;
    
    %.. Reduced Transition Matrix
    % nbar = 3;
    % aa = exp(-1/(c_lv-1));
    % Pbar0 = zeros(nbar*nx);
    % Pbar0(1:nx,1:nx) = Cp*Pq1*Pf;
    % Pbar0(1:nx, end-(nx-1):end) = rho*Cp*Pq1*Pq2*Pf;
    % 
    % Pbar0(nx+1:2*nx,1:nx) = Cm*Pq1*Pf;
    % Pbar0(nx+1:2*nx,nx+1:2*nx) = aa*Pq1;
    % Pbar0(nx+1:2*nx, end-(nx-1):end) = rho*Cm*Pq1*Pq2*Pf;
    % 
    % Pbar0(2*nx+1:3*nx, 2*nx+1:3*nx) = (1-aa)*Pq1*Pf^(c_lv-1);
    % Pbar0(2*nx+1:3*nx, 2*nx+1:3*nx) = (1-rho)*Pq1*Pf;
    % 
    % Pbar1 = zeros(nbar*nx);
    % Pbar1(1:nx,1:nx) = Cp*Pf;
    % Pbar1(1:nx, end-(nx-1):end) = rho*Cp*Pq2*Pf;
    % 
    % Pbar1(nx+1:2*nx,1:nx) = Cm*Pf;
    % Pbar1(nx+1:2*nx,nx+1:2*nx) = aa*I;
    % Pbar1(nx+1:2*nx, end-(nx-1):end) = rho*Cm*Pq2*Pf;
    % 
    % Pbar1(2*nx+1:3*nx, 1*nx+1:2*nx) = (1-aa)*Pf^(c_lv-1);
    % Pbar1(2*nx+1:3*nx, 2*nx+1:3*nx) = (1-rho)*Pf;
    % 
    % P1 = Pbar1^(c_plane-1);
    % P00 = Pbar0*P1;
    % 
    % % Prob. Dist right after RAAN Contact
    % PI_0 = limitdist(P00'); 
    % PI_hr = zeros(nbar*nx, c_plane); % Seperate Vector
    % PI_hr(:,1) = PI_0;
    % for i = 2:c_plane
    %     PI_hr(:,i) = Pbar1*PI_hr(:,i-1);
    % end

    ntp = (c_lv - 1);
    dtp = floor(ntp/d_cut);
    ctp = ntp - dtp*d_cut;
    mm = dtp*ones(1,d_cut);
    mm(1:ctp) = mm(1:ctp) + 1;
    nbar = 2 + d_cut;
    Pbar0 = zeros(nbar*nx);
    Pbar1 = zeros(nbar*nx);

    for i = 1:nbar
        if i == 1
            Pbar0(1:nx,1:nx) = Cp*Pq1*Pf;
            Pbar0(1:nx, end-(nx-1):end) = rho*Cp*Pq1*Pq2*Pf;

            Pbar1(1:nx,1:nx) = Cp*Pf;
            Pbar1(1:nx, end-(nx-1):end) = rho*Cp*Pq2*Pf;
        elseif i == 2
            Pbar0(nx+1:2*nx,1:nx) = Cm*Pq1*Pf;
            Pbar0(nx+1:2*nx,nx+1:2*nx) = (mm(1)-1)/mm(1)*Pq1;
            Pbar0(nx+1:2*nx, end-(nx-1):end) = rho*Cm*Pq1*Pq2*Pf;

            Pbar1(nx+1:2*nx,1:nx) = Cm*Pf;
            Pbar1(nx+1:2*nx,nx+1:2*nx) = (mm(1)-1)/mm(1)*I;
            Pbar1(nx+1:2*nx, end-(nx-1):end) = rho*Cm*Pq2*Pf;
        elseif i == nbar
            Pbar0((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = 1/mm(i-2)*Pq1*Pf^(mm(i-2));
            Pbar0(end-(nx-1):end, end-(nx-1):end) = (1-rho)*Pq1*Pf;

            Pbar1((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = 1/mm(i-2)*Pf^(mm(i-2));
            Pbar1(end-(nx-1):end, end-(nx-1):end) = (1-rho)*Pf;
        else
            Pbar0((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = 1/mm(i-2)*Pq1*Pf^(mm(i-2));
            Pbar0((i-1)*nx+1:i*nx, (i-1)*nx+1:i*nx) = (mm(i-1)-1)/mm(i-1)*Pq1;

            Pbar1((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = 1/mm(i-2)*Pf^(mm(i-2));
            Pbar1((i-1)*nx+1:i*nx, (i-1)*nx+1:i*nx) = (mm(i-1)-1)/mm(i-1)*I;
        end
    end
    P1 = Pbar1^(c_plane-1);
    P00 = Pbar0*P1;

    % Prob. Dist right after RAAN Contact
    PI_0 = limitdist(P00'); 
    PI_hr = zeros(nbar*nx, c_plane); % Seperate Vector
    PI_hr(:,1) = PI_0;
    for i = 2:c_plane
        PI_hr(:,i) = Pbar1*PI_hr(:,i-1);
    end

    % ntp = (c_lv - 1);
    % dtp = floor(ntp/d_cut);
    % ctp = ntp - dtp*d_cut;
    % mm = dtp*ones(1,d_cut);
    % mm(1:ctp) = mm(1:ctp) + 1;
    % nbar = 2 + d_cut;
    % Pbar0 = zeros(nbar*nx);
    % Pbar1 = zeros(nbar*nx);
    % 
    % for i = 1:nbar
    %     if i == 1
    %         Pbar0(1:nx,1:nx) = Cp*Pq1*Pf;
    %         Pbar0(1:nx, end-(nx-1):end) = rho*Cp*Pq1*Pq2*Pf;
    % 
    %         Pbar1(1:nx,1:nx) = Cp*Pf;
    %         Pbar1(1:nx, end-(nx-1):end) = rho*Cp*Pq2*Pf;
    %     elseif i == 2
    %         Pbar0(nx+1:2*nx,1:nx) = Cm*Pq1*Pf;
    %         Pbar0(nx+1:2*nx, end-(nx-1):end) = rho*Cm*Pq1*Pq2*Pf;
    % 
    %         Pbar1(nx+1:2*nx,1:nx) = Cm*Pf;
    %         Pbar1(nx+1:2*nx, end-(nx-1):end) = rho*Cm*Pq2*Pf;
    %     else
    %         Pbar0((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = Pq1*Pf^(mm(i-2));
    %         Pbar1((i-1)*nx+1:i*nx, (i-2)*nx+1:(i-1)*nx) = Pf^(mm(i-2));
    %     end
    % end
    % Pbar0(end-(nx-1):end, end-(nx-1):end) = (1-rho)*Pq1*Pf;
    % Pbar1(end-(nx-1):end, end-(nx-1):end) = (1-rho)*Pf;
    % P1 = Pbar1^(c_plane-1);
    % P00 = Pbar0*P1;
    % 
    % % Prob. Dist right after RAAN Contact
    % PI_0 = limitdist(P00'); 
    % PI_hr = zeros(nbar*nx, c_plane); % Seperate Vector
    % PI_hr(:,1) = PI_0;
    % for i = 2:c_plane
    %     PI_hr(:,i) = Pbar1*PI_hr(:,i-1);
    % end
    
    % Reconstruct Full Distribution
    % nbar_f = 2 + (c_lv - 1);
    % PI_hrf = zeros(nbar_f*nx, c_plane); % Seperate Vector
    % PI_hrf(1:2*nx,:) = PI_hr(1:2*nx,:);
    % PI_hrf(end-(nx-1):end,:) = PI_hr(end-(nx-1):end,:);
    % 
    % for i = 3:(c_lv - 1)
    %     idx_k = (i-2)*nx+1:(i-1)*nx;  
    %     idx_kk = (i-1)*nx+1:i*nx;
    %     for j = 1:c_plane
    %         if j == 1
    %             PI_hrf(idx_kk,j) = Pq1*Pf*PI_hrf(idx_k,end);
    %         else
    %             PI_hrf(idx_kk,j) = Pf*PI_hrf(idx_k,j-1);
    %         end
    %     end
    % end
    % 
    % %.. Normalize
    % for j = 1:c_plane
    %     PI_hrf(:,j) = PI_hrf(:,j)/sum(PI_hrf(:,j));
    % end
    
    % Prob. Dist during the RAAN Contact Cycle
    pi_hr = zeros(nx, c_plane); % Summed Vector
    for i = 1:c_plane
        pi_hr(:,i) = sum(reshape(PI_hr(:,i),nx,[]),2);
    end
    

    % Prob. Dist right before RAAN Contact
    pi_r1 = Pf*pi_hr(:,end);
    
    % Prob. Dist right after RAAN Contact
    pi_q1 = Pq1*pi_r1;
    
    % Demand Distribution
    pi_dmd = zeros(Dmax+1,1);
    pi_dmd(1) = sum(pi_r1(1:find(x==R1+1))); % Prob of zero order: X = R1+1 ~ Xmax
    for i = 2:Dmax+1
        % Prob of (i-1) order: X = R1+1-Q1*(i-1) ~ R1-Q1*(i-2)
        % ex) i=2 -> 1 order: X = R1+1-Q1 ~ R1
        x_ref0 = R1-Q1*(i-2); % Start Idx for (i-1) order
        x_reff = R1+1-Q1*(i-1); % End Idx for (i-1) order
        if x_ref0 < 0 % If idx out of range
            pi_dmd = pi_dmd(1:i-1);
            break
        elseif x_reff >= 0
            pi_dmd(i) = sum(pi_r1(find(x==x_ref0):find(x==x_reff)));
        elseif x_reff < 0
            pi_dmd(i) = sum(pi_r1(find(x==x_ref0):end));
        end
    end

    % Prob. Dsit for Direct Resupply
    PI_Q2 = PI_hr(end-nx+1:end,:); % pi_wp
    PI_R2 = PI_hr(nx+1:2*nx,:); % pi_tp1
    for i = 1:c_plane
        PI_Q2(:,i) = rho*Pq2*Pf*PI_Q2(:,i); %pi_q2 = rho*Pq2*Pf*pi_wp
    end
    pi_q2 = sum(PI_Q2,2);
    pi_q2 = pi_q2/sum(pi_q2);
    pi_r2 = sum(PI_R2,2);
    pi_r2 = pi_r2/sum(pi_r2);
    
    for i = 1:c_plane
        PI_Q2(:,i) = PI_Q2(:,i)/sum(PI_Q2(:,i));
        PI_R2(:,i) = PI_R2(:,i)/sum(PI_R2(:,i));
    end
    
    %.. Output
    PI.PI_hr = PI_hr;
    PI.pi_hr = pi_hr;
    
    PI.pi_q1 = pi_q1;
    PI.pi_r1 = pi_r1;
    PI.pi_dmd = pi_dmd;
    
    PI.PI_Q2 = PI_Q2;
    PI.PI_R2 = PI_R2;
    PI.pi_q2 = pi_q2;
    PI.pi_r2 = pi_r2;
end


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