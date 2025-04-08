function [J, Cost] = CostHybridResupply(X, ParaCost, ParaInPlane, ParaParking)
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
    
    %.. Analyze the performance
    [PI_i, PI_p, T_i, T_p] = SolveHybridProb(100,ParaInPlane,ParaParking);
    
    %.. Cost per unit time step for each different list
    C_build = ManufacturingCost(T_i.T_q2, T_p.T_avg, ParaCost, ParaInPlane, ParaParking);
    C_hold = HoldingCost(PI_i.X, PI_i.pi_avg, PI_p.X, PI_p.pi_avg, ParaCost, ParaInPlane, ParaParking);
    C_launch = LaunchCost(T_i.T_q2, T_p.T_avg, ParaCost, ParaInPlane, ParaParking);
    C_transfer = TransferCost(PI_i.X, PI_i.pi_r1, PI_i.pi_q1, T_i.T_q1, ParaCost, ParaInPlane, ParaParking);
    
    %.. Total Cost
    J = C_build + C_hold + C_launch + C_transfer;
    
    %.. Save cost for each term
    Cost.C_build = C_build;
    Cost.C_hold = C_hold;
    Cost.C_launch = C_launch;
    Cost.C_transfer = C_transfer;
end

%% Manufacturing Cost function
function [C_build] = ManufacturingCost(Ti_hr2, Tp_hr, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    c_build = ParaCost.c_build;
    N_plane = ParaInPlane.N_orbit;
    N_park = ParaParking.N_orbit;
    qi1 = ParaInPlane.Q1;
    qi2 = ParaInPlane.Q2;
    qp = ParaParking.Q;
    
    % Manufacturing Cost per unit time step
    C_build = c_build*( (qi2*N_plane)/Ti_hr2 + (qi1*qp*N_park)/Tp_hr);
end

%% Holding Cost function
function [C_hold] = HoldingCost(Xi, P_Xi, Xp, P_Xp, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    c_hold_plane = ParaCost.c_hold_plane;
    c_hold_park = ParaCost.c_hold_park;
    N_plane = ParaInPlane.N_orbit;
    N_park = ParaParking.N_orbit;
    N_sat = ParaInPlane.N_sat;
    qi1 = ParaInPlane.Q1;
    
    % Avg. # of Spares in in-plane orbits
    Si_k = Xi - N_sat; % (Xi_max, ... , 1, 0) - N_sat
    Si_k(Si_k <= 0) = 0;
    Si_avg = dot(Si_k,P_Xi);
    
    % Avg. # of Spares in parking orbits
    Sp_avg = qi1*dot(Xp,P_Xp);
    
    % Holding Cost per unit time step
    C_hold = c_hold_plane*N_plane*Si_avg + c_hold_park*N_park*Sp_avg;
end

%% Launch Cost function
function [C_launch] = LaunchCost(Ti_hr2, Tp_hr, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    c_lv_small_full = ParaCost.c_lv_small_full;
    c_lv_small_part = ParaCost.c_lv_small_part;
    c_lv_heavy_full = ParaCost.c_lv_heavy_full;
    c_lv_heavy_part = ParaCost.c_lv_heavy_part;
    N_plane = ParaInPlane.N_orbit;
    N_park = ParaParking.N_orbit;
    qi1 = ParaInPlane.Q1;
    qi2 = ParaInPlane.Q2;
    qp = ParaParking.Q;
    
    % Apply minimum cost among full / partial cost
    c_launch_small = min(c_lv_small_full, qi2*c_lv_small_part);
    c_launch_heavy = min(c_lv_heavy_full, qi1*qp*c_lv_heavy_part);
    
    % Launch Cost per unit time step
    C_launch = (N_park/Tp_hr)*c_launch_heavy + (N_plane/Ti_hr2)*c_launch_small;
end

%%
function [C_trn] = TransferCost(Xi, P_Xir, P_Xiq, T_plane, ParaCost, ParaInPlane, ParaParking)
    % Input definition
    h_plane = ParaInPlane.alt;
    N_plane = ParaInPlane.N_orbit;
    qi1 = ParaInPlane.Q1;
    h_park = ParaParking.alt;

    c_fuel = ParaCost.c_fuel;
    c_trn = ParaCost.c_transfer;
    m_sat = ParaCost.m_sat;
    m_bus = ParaCost.m_bus;
    Vex = ParaCost.Vex;
    
    % Fuel for both (bus + spares) transfer
    m_fuel = ComputeTransferFuelPerBatch(h_park, h_plane, m_bus, m_sat, qi1, Vex);
    
    % Avg. # of transfer for every review period (T_plane)
    N_trn = (dot(Xi,P_Xiq) - dot(Xi,P_Xir))/qi1; % [batch]
    
    % Transfer cost per unit time step
    C_trn = N_plane/T_plane * N_trn * ( c_fuel*m_fuel + c_trn );
end

function [m_fuel] = ComputeTransferFuelPerBatch(h_park, h_plane, m_bus, m_sat, q1i, Vex)
    global R_earth mu_earth
    ai = R_earth + h_plane;
    ap = R_earth + h_park;
    DV = sqrt(mu_earth/ap)*( sqrt(2*ai/(ai+ap)) - 1 ) + sqrt(mu_earth/ai)*( 1 - sqrt(2*ap/(ai+ap)) );
    m_fuel = (m_bus + q1i*m_sat)*( exp(DV/Vex) - 1 );
end
