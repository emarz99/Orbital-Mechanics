function [a, e, i, OM, om, th] = car2kep(r, v, mu)

% car2kep - Conversion from Cartesian coordinates to Keplerian elements

% INPUT:
% r    [3x1] Position vector [km]
% v    [3x1] Velocity vector [km/s]
% mu   [1x1] Gravitational parameter [km^3/s^2]

% OUTPUT:
% a    [1x1] Semi-major axis [km]
% e    [1x1] Eccentricity [-]
% i    [1x1] Inclination [rad]
% OM   [1x1] RAAN [rad]
% om   [1x1] Pericentre anomaly [rad]
% th   [1x1] True anomaly [rad]

r_m = norm(r);   % Absolute value of the position vector [km]
v_m = norm(v);   % Absolute value of the velocity vector [km/s]

h = cross(r,v);  % Angular momentum [km^2/s]

i = acos( h(3,1)/ norm(h));   % Inclination [rad]

e_v = (1/mu)*[(v_m^2-mu/r_m)*r-dot(r,v)*v];   % Eccentricity vector [-]
e = norm(e_v);   % Eccentricity [-]

E = 0.5*v_m^2-(mu/r_m);   % Specific energy [km^2/s^2]
a = -mu/(2*E);   % Semimajor axis [km]

K = [0; 0; 1]; 
N = cross(K,h);   % Node line vector [-]
N_m = norm(N);    % Absolute value of Node line vector [-]

if N(2,1)>=0
    OM = acos(N(1,1)/N_m);   % Right ascension of the ascendent node [rad]
else
    OM = 2*pi-acos(N(1,1)/N_m);
end
if e_v(3,1)>=0
    om = acos((dot(N,e_v))/(N_m*e));   % Argument of perigree [rad]
else
    om = 2*pi-acos(dot(N,e_v)/(N_m*e));
end

v_r = dot(r,v)/r_m;   % Radial velocity [km/s]

if v_r>=0
    th = acos(dot(e_v,r)/(e*r_m));   % True anomaly [rad]
else
    th = 2*pi-acos(dot(e_v,r)/(e*r_m));
end
end








