% Assignment 2 Data

a = 2.4039*1e+4;             % Initial semi-major axis [km]
e = 0.7236;                  % Initial eccentricity [-]
i = deg2rad(4.1093);         % Initial inclination of the orbital plane [rad]
OM = deg2rad(30);            % Initial right ascension of the ascending node [rad]
om = deg2rad(60);            % Initial argument of pericentre [rad]
th0 = pi/2;                  % Initial true anomaly [rad]
h_p = 273.370;               % Initial pericentre altitude [km]
ratio = 7/3;                 % Ratio between the number of revolution of the satellite
                             % and the number of rotations of the planet [-]
cD = 2.1;                    % Drag coefficient [-]
A_M = 0.0054/(1e6);          % Cross area to mass ratio (for the drag computation) [km^2/kg]
mu_E = astroConstants(13);   % Planetary constant of the Earth [km^3/s^2]
R = astroConstants(23);      % Mean radius of the Earth
J2 = astroConstants(9);      % Zonal harmonics (two dimensions) [-]
om_E = deg2rad(15.04)/3600;  % [rad/s] Earth's rotation velocity
