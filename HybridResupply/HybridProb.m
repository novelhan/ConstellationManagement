% Def
% Compute Exact State Distrubution of Hybrid Resupply Method
%
% Input
% iter_max: Max iteration number for fixed point iteration method
% ParaInPlane: In-Plane orbit parameter structure
%   ParaFail: Structure for Markov model
%   ParaDirect: Structure for Direct channel 
%   ParaIndirect: Structure for Indirect channel 
%   ParaDim: Structure for Dim. reduction
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


function [PI_i, PI_p, T_i, T_p, err] = HybridProb(iter_max, ParaInPlane, ParaParking)
    %.. Assume Perfect Availability for initialization
    [~, PI_i] = HybridInPlane(ParaInPlane.ParaFail, ParaInPlane.ParaDirect,...
                              ParaInPlane.ParaIndirect, ParaInPlane.ParaDim);
    if ParaParking.method == 0 % Select solution method for parking orbit analysis
    [~, PI_p] = ExactInDirectParkTime(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                      ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
    else
    [~, PI_p] = ExactInDirectParkRatio(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                       ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
    end
    
    %.. Fixed Point Root Finding
    xx_pre = [mean(PI_i.pi_hr,2); PI_p.pi_ir];
    err = zeros(1,iter_max);
    for iter = 1:iter_max
        %.. Update the solution
        ParaInPlane.ParaIndirect.kappa = PI_p.Pdav;
        [~, PI_i, T_i] = HybridInPlane(ParaInPlane.ParaFail, ParaInPlane.ParaDirect,...
                                       ParaInPlane.ParaIndirect, ParaInPlane.ParaDim);
        if ParaParking.method == 0 % Select solution method for parking orbit analysis
        [~, PI_p, T_p] = ExactInDirectParkTime(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                               ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
        else
        [~, PI_p, T_p] = ExactInDirectParkRatio(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                                ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
        end
        
        %.. Check Convergence
        xx = [mean(PI_i.pi_hr,2); PI_p.pi_ir];
        err(iter) = norm(xx - xx_pre);
%         disp(['Iteration:', num2str(iter), ',  Error:', num2str(err(iter))])
        if err(iter) < 10^(-7)
            err = err(1:iter);
            break;
        else
            xx_pre = xx;
        end
    end
end