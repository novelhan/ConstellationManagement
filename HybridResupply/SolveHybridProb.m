% Def
% Compute Exact State Distrubution of Hybrid Resupply Method
%
% Input
% iter_max: Max iteration number for fixed point iteration method
% ParaInPlane: In-Plane orbit parameter structure
% ParaParking: Parking orbit parameter structure
% The detailed information of each structure can be found in each sub-function
%
% Output
% PI_i: Set of Stationary State Distribution for in-plane orbit
% PI_p: Set of Stationary State Distribution for parking orbit
% T_i: Set of time duration for cycles for in-plane orbit
% T_p: Set of time duration for cycles for parking orbit
% err: residual history of fixed point iteration
%
% Reference


function [PI_i, PI_p, T_i, T_p, err] = SolveHybridProb(iter_max, ParaInPlane, ParaParking)
    %.. Assume Perfect Availability for initialization
    nx_i = ParaInPlane.Q1 + ParaInPlane.Q2 + max(ParaInPlane.R1, ParaInPlane.R2) + 1;
    Pdav = 0.9*ones(1,nx_i);
    [PI_i] = SolveHybridInPlane(Pdav, ParaInPlane);
	[PI_p] = SolveInDirectPark(PI_i.pi_dmd, ParaParking);
    
    %.. Fixed Point Root Finding
    xx_pre = [PI_i.pi_avg; PI_p.pi_avg];
    err = zeros(1,iter_max);
    for iter = 1:iter_max
        %.. Update the solution
        [PI_i, T_i] = SolveHybridInPlane(PI_p.Pdav, ParaInPlane);
        [PI_p, T_p] = SolveInDirectPark(PI_i.pi_dmd, ParaParking);
        
        %.. Check Convergence
        xx = [PI_i.pi_avg; PI_p.pi_avg];
        err(iter) = norm(xx - xx_pre);
%         disp(['Iteration:', num2str(iter), ',  Error:', num2str(err(iter))])
        if err(iter) < 10^(-6)
            err = err(1:iter);
            break;
        else
            xx_pre = xx;
        end
    end
end