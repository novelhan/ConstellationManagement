function [c,ceq] = ConstInDirectResupply(X, ParaConst, ParaInPlane, ParaParking)
    %.. Set Design Variable
    ParaInPlane.Q = X(1); % Qi
    ParaInPlane.R = X(2); % Ri
    ParaParking.Q = X(3); % Qp
    ParaParking.R = X(4); % Rp
    ParaParking.N_orbit = X(5); % N_park
    ParaParking.alt = X(6); 
    
    %.. Compute RAAN Drift Period
    [dt_plane, dt_park] = ComputeRaanPeriod(ParaInPlane, ParaParking);
    ParaInPlane.dt_plane = dt_plane;
    ParaParking.dt_park = dt_park;

    %.. Analyze the performance
    [PI_i, PI_p] = SolveInDirectProb(100,ParaInPlane,ParaParking);
    
    %.. InPlane Orbit Resilience Constraint: P(X < Xref) <= eps
    Si_k = PI_i.X - ParaInPlane.N_sat; % (Xi_max, ... , 1, 0) - N_sat
    c(1) = sum(PI_i.pi_avg(Si_k<0)) - ParaConst.p_loss_plane; 
    
    %.. Parking Orbit Resilience Constraint: P(X == 0) <= eps
    c(2) = sum(PI_p.pi_avg(end)) - ParaConst.p_loss_park; 
    
    %.. Constraint on max # of spares for single LV: Qi*Qp <= Q_max
    c(3) = X(1) * X(3) - ParaConst.N_lv_max_heavy;
    
    %.. For time based analysis: Rp < Qp => Rp + 1 <= Qp
    if ParaParking.method == 0
        c(4) = X(4) + 1 - X(3);
    end
    
    %.. Lower and Upper bound for Design Variable
    % Handled by external feature
    
    %.. No equality constraint
    ceq = [];
end