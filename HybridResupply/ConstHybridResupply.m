function [c,ceq] = ConstHybridResupply(X, ParaConst, ParaInPlane, ParaParking)
    %.. Set Design Variable
    ParaInPlane.Q1 = X(1); % Qi1
    ParaInPlane.R1 = X(2); % Ri1
    ParaInPlane.Q2 = X(3); % Qi2
    ParaInPlane.R2 = X(4); % Ri2
    ParaParking.Q = X(5); % Qp
    ParaParking.R = X(6); % Rp
    ParaParking.N_orbit = X(7); % N_park
    ParaParking.alt = X(8); 
    
    %.. Compute RAAN Drift Period
    [dt_plane, dt_park] = ComputeRaanPeriod(ParaInPlane, ParaParking);
    ParaInPlane.dt_plane = dt_plane;
    ParaParking.dt_park = dt_park;

    
    if  ParaParking.method == 0 && X(6) + 1 - X(5) > 0
        c(1:5) = 1;
    else
        
    %.. Analyze the performance
    [PI_i, PI_p] = SolveHybridProb(100,ParaInPlane,ParaParking);
    
    %.. InPlane Orbit Resilience Constraint: P(X < Xref) <= eps
    Si_k = PI_i.X - ParaInPlane.N_sat; % (Xi_max, ... , 1, 0) - N_sat
    c(1) = sum(PI_i.pi_avg(Si_k<0)) - ParaConst.p_loss_plane; 
    
    %.. Parking Orbit Resilience Constraint: P(X == 0) <= eps
    c(2) = sum(PI_p.pi_avg(end)) - ParaConst.p_loss_park; 
    
    %.. Constraint on max # of spares for single heavy LV: Qi1*Qp <= Q1_max
    c(3) = X(1) * X(5) - ParaConst.N_lv_max_heavy;
    
    %.. Constraint on max # of spares for single small LV: Qi2 <= Q2_max
    c(4) = X(3) - ParaConst.N_lv_max_small;
    
    %.. For time based analysis: Rp < Qp => Rp + 1 <= Qp
    if ParaParking.method == 0
        c(5) = X(6) + 1 - X(5);
    end
    end
    %.. Lower and Upper bound for Design Variable
    % Handled by external feature
    
    %.. No equality constraint
    ceq = [];
end