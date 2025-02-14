%% RAAN Period Computing Function
function [dt_plane, dt_park] = ComputeRaanPeriod(ParaInPlane, ParaParking)
    %.. Earth Parameters
    global R_earth mu_earth J2_earth
    
    %.. In-Orbit Parameters
    a_plane     =   R_earth + ParaInPlane.alt;
    e_plane     =   0;
    i_plane     =   ParaInPlane.inc;
    n_plane     =   sqrt(mu_earth/a_plane^3);
    N_plane     =   ParaInPlane.N_orbit;
    Wdot_orbit  =   -3/2 * J2_earth* n_plane * R_earth^2 / (a_plane*(1-e_plane))^2 * cos(i_plane);
    
    %.. Parking-Orbit Parameters
    a_park      =   R_earth + ParaParking.alt;
    e_park      =   0;
    i_park      =   ParaParking.inc;
    n_park      =   sqrt(mu_earth/a_park^3);
    N_park      =   ParaParking.N_orbit;
    Wdot_park   =   -3/2 * J2_earth * n_park * R_earth^2 / (a_park*(1-e_park))^2 * cos(i_park);

    %.. RAAN Drift period for in-plane and parking orbit
    Wdrift      =   abs(Wdot_park - Wdot_orbit);
    dt_plane    =   2*pi/ N_park / Wdrift / 24 / 3600; % RAAN Lead Time [day]
    dt_park     =   2*pi/ N_plane / Wdrift / 24 / 3600; % RAAN Lead Time [day]
end