% Def
% Compute Exact State Distrubution of Indirect Resupply Method
%
% Input
% iter_max: Max iteration number for fixed point iteration method
% ParaInPlane: In-Plane orbit parameter structure
% ParaParking: Parking orbit parameter structure
% The detailed information of each structure can be found in 
% ExactInDirectPlane.m and ExactInDirectPark.m
%
% Output
% PI_i: Set of Stationary State Distribution for in-plane orbit
% PI_p: Set of Stationary State Distribution for parking orbit
% T_i: Set of time duration for cycles for in-plane orbit
% T_p: Set of time duration for cycles for parking orbit
% err: residual history of fixed point iteration
%
% Reference
% Analysis and Design of Satellite Constellation Spare Strategy Using Markov Chain
% https://doi.org/10.48550/arXiv.2408.09250


function [PI_i, PI_p, T_i, T_p, err] = ExactInDirectProb(iter_max, ParaInPlane, ParaParking)
    %.. Assume Perfect Availability for initialization
    nx_i = ParaInPlane.Q + ParaInPlane.R + 1;
    Pdav = 0.9*ones(1,nx_i);
    [~, PI_i] = ExactInDirectPlane(ParaInPlane.f_sim, ParaInPlane.f_type, ParaInPlane.n_sat, Pdav,...
                                   ParaInPlane.Q, ParaInPlane.R, ParaInPlane.dt_mc, ParaInPlane.dt_plane);
    if ParaParking.method == 0 % Select solution method for parking orbit analysis
    [~, PI_p] = ExactInDirectParkTime(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                      ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
    else
    [~, PI_p] = ExactInDirectParkRatio(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                       ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
    end
    
    %.. Fixed Point Root Finding
    xx_pre = [PI_i.pi_ir; PI_p.pi_ir];
    err = zeros(1,iter_max);
    for iter = 1:iter_max
        %.. Update the solution
        [~, PI_i, T_i] = ExactInDirectPlane(ParaInPlane.f_sim, ParaInPlane.f_type, ParaInPlane.n_sat, PI_p.Pdav,...
                                            ParaInPlane.Q, ParaInPlane.R, ParaInPlane.dt_mc, ParaInPlane.dt_plane);
        if ParaParking.method == 0 % Select solution method for parking orbit analysis
        [~, PI_p, T_p] = ExactInDirectParkTime(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                               ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
        else
        [~, PI_p, T_p] = ExactInDirectParkRatio(PI_i.pi_dmd, ParaParking.Q, ParaParking.R, ParaParking.dt_mc, ...
                                                ParaParking.dt_park, ParaParking.mu_lv, ParaParking.dt_lv);
        end
        
        %.. Check Convergence
        xx = [PI_i.pi_ir; PI_p.pi_ir];
        err(iter) = norm(xx - xx_pre);
%         disp(['Iteration:', num2str(iter), ',  Error:', num2str(err(iter))])
        if err(iter) < 10^(-5)
            err = err(1:iter);
            break;
        else
            xx_pre = xx;
        end
    end
end