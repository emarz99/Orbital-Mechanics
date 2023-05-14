function [r,v] = kep2car(a, e, i, OM, om, th, mu)

% kep2car - Conversion from Keplerian elements to Cartesian coordinates

% INPUT:
% a    [1x1] Semi-major axis [km]
% e    [1x1] Eccentricity [-]
% i    [1x1] Inclination [rad]
% OM   [1x1] RAAN [rad]
% om   [1x1] Pericentre anomaly [rad]
% th   [1x1] True anomaly [rad]
% mu   [1x1] Gravitational parameter [km^3/s^2]

% OUTPUT:
% r    [3x1] Position vector [km]
% v    [3x1] Velocity vector [km/s]

p = a*(1-e^2);   % [1x1] semilatus rectum [km]
r_m = p/(1+e*cos(th));   % Absolute value of the position vector [km]
r_pf = r_m*[cos(th); sin(th); 0];   % Position vector (perifocal coordinate system) [km]
v_pf = sqrt(mu/p)*[-sin(th); e+cos(th); 0];   % Velocity vector (perifocal coordinate system) [km/s]

% Switch from the perifocal coordinate system to the Geocentric equatorial coordinates (i,j,k):

R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];   % Rotation matrix for OM around k axis
R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];   % Rotation matrix for i around i' axis [-]
R3_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];   % Rotation matrix for om around k'' axis

T = R3_om*R1_i*R3_OM;   %matrix of trasformation PF --> ECI

r = T'*r_pf;   % Position vector (Geocentric equatorial coordinates) [km]
v = T'*v_pf;   % Velocity vector (Geocentric equatorial coordinates) [km/s]
end






