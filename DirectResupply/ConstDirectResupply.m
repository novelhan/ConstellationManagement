function [c,ceq] = ConstDirectResupply(X, ParaConst, ParaInPlane)
    %.. Set Design Variable
    ParaInPlane.Q = X(1);
    ParaInPlane.R = X(2);

    %.. Analyze the performance
    [PI, ~] = SolveDirectProb(ParaInPlane);
    
    %.. Resilience Constraint: P(X < Xref) <= eps
    Si_k = PI.X - ParaInPlane.N_sat; % (Xi_max, ... , 1, 0) - N_sat
    c = sum(PI.pi_dr(Si_k<0)) - ParaConst.p_loss; 
     
    %.. Lower and Upper bound for Design Variable
    % Handled by external feature
    
    %.. No equality constraint
    ceq = [];
end