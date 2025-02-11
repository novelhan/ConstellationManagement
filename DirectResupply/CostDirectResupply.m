function [J, Cost] = CostDirectResupply(X, ParaCost, ParaInPlane)
    %.. Set Design Variable
    ParaInPlane.Q = X(1);
    ParaInPlane.R = X(2);

    %.. Analyze the performance
    [PI, T] = SolveDirectProb(ParaInPlane);
    
    %.. Cost per unit time step for each different list
    C_build = DirectManufacturingCost(T.T_dr, ParaCost, ParaInPlane);
    C_hold = DirectHoldingCost(PI.Xi, PI.pi_dr, ParaCost, ParaInPlane);
    C_launch = DirectLaunchCost(T.T_dr, ParaCost, ParaInPlane);
    
    %.. Total Cost
    J = C_build + C_hold + C_launch;
    
    %.. Save cost for each term
    Cost.C_build = C_build;
    Cost.C_hold = C_hold;
    Cost.C_launch = C_launch;
end

%% Manufacturing Cost function
function [C_build] = DirectManufacturingCost(T, ParaCost, ParaInPlane)
    % Input definition
    c_build = ParaCost.c_build;
    N_plane = ParaInPlane.N_plane;
    q = ParaInPlane.Q;
    
    % Manufacturing Cost per unit time step
    C_build = c_build*N_plane*q/T;
end

%% Holding Cost function
function [C_hold] = DirectHoldingCost(Xi, P_Xi, ParaCost, ParaInPlane)
    % Input definition
    c_hold = ParaCost.c_hold;
    N_plane = ParaInPlane.N_plane;
    N_sat = ParaInPlane.N_sat;
    
    % Avg. # of Spares
    Si_k = Xi - N_sat; % (Xi_max, ... , 1, 0) - N_sat
    Si_k(Si_k <= 0) = 0;
    Si_avg = dot(Si_k,P_Xi);
    
    % Holding Cost per unit time step
    C_hold = c_hold*N_plane*Si_avg;
end

%% Launch Cost function
function [C_launch] = DirectLaunchCost(T, ParaCost, ParaInPlane)
    % Input definition
    c_lv_small_full = ParaCost.c_lv_small_full;
    c_lv_small_part = ParaCost.c_lv_small_part;
    N_plane = ParaInPlane.N_plane;
    q = ParaInPlane.Q;
    
    % Apply minimum cost among full / partial cost
    c_launch_small = min(c_lv_small_full, q*c_lv_small_part);
    
    % Launch Cost per unit time step
    C_launch = (N_plane/T)*c_launch_small;
end


