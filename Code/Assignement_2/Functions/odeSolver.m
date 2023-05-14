function [time, kep]=odeSolver(kep0, mu, J2, R, A_M, cD, om_E, tspan, options)

% odeSolver - Solver for the ordinary differential equations in Keplerian
% elements with respect to time

% INPUT:
% kep0     Initial conditions in Keplerian elements
% mu       Planetary constants of the planet [km^3/s^2]
% R        Mean radius of the planet [km]
% J2       Zonal harmonics (two dimensions) [-]
% A_M      Cross area to mass ratio (for the drag computation) [km^2/kg]
% cD       Drag coefficient (for the drag computation) [-]
% om_E     Angular velocity of the planet along h [rad/s]
% tspan    Time intervall [s]
% options  Options for the ode solver

% OUTPUT:
% time    [1xN] time vector [s]
% kep     [Nx6] time history Keplerian elements 

% AUTHORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

[time, kep] = ode113(@(t,kep) gaussEoM_rsw( t, kep, mu, a_Per_rsw(t, kep, mu, R, J2, A_M, cD, om_E) ), tspan, kep0, options);

end