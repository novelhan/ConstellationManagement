% Def
% Compute Exact State Distrubution of InPlane Orbit Spares under Hybrid Resupply Method
%
% Input
% Failure Related Input Strcture
%   dt_mc: unit time step
%   f_ref: failure rate per day (reference)
%   f_type: failure type
%   n_sat: # of nominal satellite
%
% Direct Resupply Related Input Strcture
%   mu: exponential mean time of lead time for direct resupply
%   dt_lv: process time of lead time for direct resupply
%   Q2: Reorder quantity for direct resupply
%   R2: Reorder level for direct resupply
%
% Indirect Resupply Related Input Strcture
%   dt_plane: review period (duration unitl meet subsequent parking orbit)
%   kappa: Parking orbit spare availability probability distribution
%   Q1: Reorder quantity for indirect resupply
%   R1: Reorder level for indirect resupply
%
% Dimension Reduction Related Input Strcture
%   dim_flag: flag to decide using dimension reduction (0:No reduction, 1: Apply reduction)
%   x_min: The number of state dimension to cut
%   n_seg: The number of segments for approximating LV processing period
%
% Output
% PI: Set of Stationary State Distribution
% T: Set of time duration for cycles
%
% Reference
%

function [PI, T] = SolveHybridInPlane(kappa, ParaInPlane)
    %.. Failure Parameter
    dt_mc = ParaInPlane.dt_mc;
    f_mc = ParaInPlane.f_ref * dt_mc;
    f_type = ParaInPlane.f_type;
    n_sat = ParaInPlane.N_sat;
    
    %.. Direct Resupply Parameter
    mu = ParaInPlane.mu;
    dt_lv = ParaInPlane.dt_lv;
    Q2 = ParaInPlane.Q2;
    R2 = ParaInPlane.R2;
    
    %.. Indirect Resupply Parameter
    dt_plane = ParaInPlane.dt_plane;
    Q1 = ParaInPlane.Q1;
    R1 = ParaInPlane.R1;
    
    %.. LV Param
    lam = 1/mu;
    rho = (1 - exp(-lam*dt_mc));
    c_lv = round(dt_lv/dt_mc);
    
    %.. Constant reorder period
    c_plane = round(dt_plane/dt_mc);
    
    %.. Dimension Reduction Parameter
    if ParaInPlane.dim_flag == 1
        x_min = min(R1,R2) - ceil((3*ParaInPlane.f_ref*n_sat)*(dt_lv + 2*mu));
        n_seg = round( (dt_lv/dt_mc - 1)/3 );
        x_min = max(x_min,0); % Minimum stock level
        n_seg = max( min(n_seg,c_lv-1), 1); % 1 <= # of segment <= c_lv-1
    else
        x_min = 0;
        n_seg = c_lv-1;
    end
    
    %.. State
    x_max = Q1+Q2+max(R1,R2); % Max can happen when spares from both channels are ordered
    x = (x_max:-1:x_min)';
    nx = length(x);
    
    %.. Failure Transition Matrix
    I = eye(nx);
    Pf = zeros(nx);
    for i = 1:nx
        if f_type == 0 % Constant Failure Rate
            Pf(i:end,i) = CustomPoisPdf(0:nx-i, n_sat*f_mc)';        
        else % Stock Level Dependant Failure Rate
            if (x_max+1-i) > n_sat
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, n_sat*f_mc)';        
            else
                Pf(i:end,i) = CustomPoisPdf(0:nx-i, (x_max+1-i)*f_mc)';  
            end
        end
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end
    
    %.. InDirect Resupply Transition Matrix
    nx1 = Q1+R1+1;
    mx1 = ceil(nx1/Q1);
    Pq1_m = zeros(mx1*Q1);
    
    % Handle the corner case where Qp+Rp < mx1
    if length(kappa) < mx1 
        kappa = [kappa; zeros(mx1,1)];
    end
    
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
     
    r = x_max - Q1 - R1;
    Pq1 = blkdiag(eye(r),Pq1_m); % same as [eye(r), zeros(r,mx1*Q1); zeros(mx1*Q1,r), Pq1_m];
    Pq1 = Pq1(1:nx,1:nx);
    
    %.. Direct Resupply Transition Matrix 
    Pq2 = zeros(nx);
    Pq2(1, 1:Q2) = ones(1,Q2); % Due to Cm operation, same as Pq2(1:Q2, 1:Q2) = eye(Q2);
    Pq2(1:nx-Q2, Q2+1:nx) = eye(nx-Q2);
    
    %.. Selection Matrix for (Q2,R2) policy
    Cm = blkdiag(zeros(x_max-R2),eye(R2+1));
    Cm = Cm(1:nx,1:nx);
    Cp = I - Cm;
    
    %.. Compute step size for each segment (distribute c_lv - 1 to n_seg)
    % ex. Given n_tp = 5 and n_seg = 3 -> nn_tp = [2 2 1]
    % ex. Given n_tp = 5 and n_seg = 4 -> nn_tp = [2 1 1 1]
    n_tp = c_lv - 1; % # of original full segments
    d_tp = floor( n_tp / n_seg); % min # of step for each n_seg
    r_tp = n_tp - d_tp*n_seg; % remaining segments to be allocated
    nn_tp = d_tp*ones(1,n_seg); % # of step for each segmets
    nn_tp(1:r_tp) = nn_tp(1:r_tp) + 1;
    
    
    %.. Full Transition Matrix
    %.. It is assumed that Failure -> Check Q2 arrival -> Apply (Q1,R1) -> Check R2
    %.. Order of state vector: [pi_np, pi_tp(1), pi_tp(2), ...., pi_tp(n_seg), pi_wp]
    L = 2 + n_seg;
    Pbar0 = zeros(L*nx); % Transition at Contact moment
    Pbar1 = zeros(L*nx); % Transition during non-contact moment
    
    for i = 1:L
        if i == 1
            Pbar0(1:nx,         1:1*nx) = Cp*Pq1*Pf;            % pi_np -> pi_np
            Pbar0(1:nx,(L-1)*nx+1:L*nx) = rho*Cp*Pq1*Pq2*Pf;    % pi_wp -> pi_np
            
            Pbar1(1:nx,         1:1*nx) = Cp*Pf;            % pi_np -> pi_np
            Pbar1(1:nx,(L-1)*nx+1:L*nx) = rho*Cp*Pq2*Pf;    % pi_wp -> pi_np
        elseif i == 2
            Pbar0(nx+1:2*nx,         1:1*nx) = Cm*Pq1*Pf;                   % pi_np -> pi_tp(1)
            Pbar0(nx+1:2*nx,      nx+1:2*nx) = (nn_tp(1)-1)/nn_tp(1)*Pq1;   % pi_tp(1) -> pi_tp(1)
            Pbar0(nx+1:2*nx,(L-1)*nx+1:L*nx) = rho*Cm*Pq1*Pq2*Pf;           % pi_wp -> pi_tp(1)
            
            Pbar1(nx+1:2*nx,         1:1*nx) = Cm*Pf;                   % pi_np -> pi_tp(1)
            Pbar1(nx+1:2*nx,      nx+1:2*nx) = (nn_tp(1)-1)/nn_tp(1)*I; % pi_tp(1) -> pi_tp(1)
            Pbar1(nx+1:2*nx,(L-1)*nx+1:L*nx) = rho*Cm*Pq2*Pf;           % pi_wp -> pi_tp(1)
        elseif i == L
            Pbar0((L-1)*nx+1:i*nx,(L-2)*nx+1:(L-1)*nx) = 1/nn_tp(L-2)*Pq1*Pf^(nn_tp(L-2));  % pi_tp(end) -> pi_wp
            Pbar0((L-1)*nx+1:L*nx,(L-1)*nx+1:    L*nx) = (1-rho)*Pq1*Pf;                    % pi_wp -> pi_wp

            Pbar1((L-1)*nx+1:i*nx,(L-2)*nx+1:(L-1)*nx) = 1/nn_tp(L-2)*Pf^(nn_tp(L-2));  % pi_tp(end) -> pi_wp
            Pbar1((L-1)*nx+1:L*nx,(L-1)*nx+1:    L*nx) = (1-rho)*Pf;                    % pi_wp -> pi_wp
        else
            Pbar0((i-1)*nx+1:i*nx,(i-2)*nx+1:(i-1)*nx) = 1/nn_tp(i-2)*Pq1*Pf^(nn_tp(i-2));  % pi_tp(k-1) -> pi_tp(k)
            Pbar0((i-1)*nx+1:i*nx,(i-1)*nx+1:    i*nx) = (nn_tp(i-1)-1)/nn_tp(i-1)*Pq1;     % pi_tp(k) -> pi_tp(k)
            
            Pbar1((i-1)*nx+1:i*nx,(i-2)*nx+1:(i-1)*nx) = 1/nn_tp(i-2)*Pf^(nn_tp(i-2));  % pi_tp(k-1) -> pi_tp(k)
            Pbar1((i-1)*nx+1:i*nx,(i-1)*nx+1:    i*nx) = (nn_tp(i-1)-1)/nn_tp(i-1)*I;   % pi_tp(k) -> pi_tp(k)
        end
    end
    Ptr = Pbar0*Pbar1^(c_plane-1);
    
    % Prob. Dist during the RAAN Contact Cycle
    PI_hr = zeros(L*nx, c_plane); % Seperated full state vector over time
    pi_hr = zeros(nx, c_plane); % Summed vector over time
    PI_hr(:,1) = limitdist(Ptr); % Prob. Dist right after RAAN Contact
    pi_hr(:,1) = sum(reshape(PI_hr(:,1),nx,[]),2);
    for i = 2:c_plane
        PI_hr(:,i) = Pbar1*PI_hr(:,i-1);
        pi_hr(:,i) = sum(reshape(PI_hr(:,i),nx,[]),2);
    end
    
    % Prob. Dist right before RAAN Contact
    pi_r1 = Pf*pi_hr(:,end);
    
    % Prob. Dist right after RAAN Contact
    pi_q1 = Pq1*pi_r1;
    
    % Demand Distribution
    Dmax = ceil(x_max/Q1);
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
    PI_q2 = zeros(nx,c_plane);
    PI_r2 = zeros(nx,c_plane);
    PI_tp = zeros(nx,c_plane);
    PI_np = PI_hr(         1:1*nx,:); % pi_np
    PI_wp = PI_hr((L-1)*nx+1:L*nx,:); % pi_wp
    
    for i = 1:c_plane
        if i == 1
            PI_q2(:,i) = rho*Pq2*Pf*PI_wp(:,end); 
            PI_r2(:,i) = Cm*Pq1*Pf*PI_np(:,end) + rho*Cm*Pq1*Pq2*Pf*PI_wp(:,end); 
        else
            PI_q2(:,i) = rho*Pq2*Pf*PI_wp(:,i-1); 
            PI_r2(:,i) = Cm*Pf*PI_np(:,i-1) + rho*Cm*Pq2*Pf*PI_wp(:,i-1); 
        end
        PI_tp(:,i) = sum(reshape(PI_hr(nx+1:(L-1)*nx,i),nx,[]),2);
    end
    pi_q2 = sum(PI_q2,2)/sum(sum(PI_q2,2));
    pi_r2 = sum(PI_r2,2)/sum(sum(PI_r2,2));
    
    %
    a_np = sum(mean(PI_np,2));
    a_wp = sum(mean(PI_wp,2));
    a_tp = sum(mean(PI_tp,2))/(c_lv-1)*c_lv; % Only k_lv - 1 terms are considered for the states
    
    %.. Reorder period
    T_q1 = c_plane*dt_mc; % Q1 review period
    T_tp = c_lv*dt_mc;
    T_np = a_np/a_tp*T_tp;
    T_wp = a_wp/a_tp*T_tp;
    T_q2 = T_np + T_tp + T_wp; % Q2 review period
    
    T.T_q1 = T_q1;
    T.T_q2 = T_q2;
    
    %.. Normalization
    for i = 1:c_plane
        PI_q2(:,i) = PI_q2(:,i)/sum(PI_q2(:,i));
        PI_r2(:,i) = PI_r2(:,i)/sum(PI_r2(:,i));
    end
    
    %.. Output
    PI.X = x;
    
    PI.PI_hr = PI_hr;
    PI.pi_hr = pi_hr;
    
    PI.pi_avg = mean(pi_hr,2);
    PI.pi_q1 = pi_q1;
    PI.pi_r1 = pi_r1;
    PI.pi_dmd = pi_dmd;
    
    PI.PI_q2 = PI_q2;
    PI.PI_r2 = PI_r2;
    PI.pi_q2 = pi_q2;
    PI.pi_r2 = pi_r2;
    
    
end
