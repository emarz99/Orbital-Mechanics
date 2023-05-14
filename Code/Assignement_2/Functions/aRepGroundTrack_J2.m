function a = aRepGroundTrack_J2(ratio,om_E,mu,J2,R,e,i,a0)

% aRepGroundTrack_J2: compute the semi-major axis for repeating ground
%                     track with J2 perturbations

% INPUT:
% ratio   [1x1] Ratio between number of revolution of the satellite and 
%               number of rotations of the planet [-]
% om_E    [1x1] Earth's rotation velocity [rad/s]
% mu      [1x1] Gravitational parameter [km^3/s^2]
% J2      [1x1] Zonal harmonics (two dimensions) [-]
% R       [1x1] Radius of the planet [km]
% e       [1x1] Eccentricity [-]
% i       [1x1] Inclination [rad]
% a0      [1x1] Initial guess for the semi-major axis [km]

% OUTPUT:
% a       [1x1] Semi-major axis of modified orbit [km]

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso
% Pasquariello Chiara

a_fun = @(a) 1/ratio - (om_E+(3*sqrt(mu)*J2*R^2/(2*(1-e^2)^2*a^(7/2)))*cos(i))/(sqrt(mu/a^3)-(3*sqrt(mu)*J2*R^2/(2*(1-e^2)^2*a^(7/2)))*((5/2)*sin(i)^2-2)+(3*sqrt(mu)*J2*R^2/(2*(1-e^2)^(3/2)*a^(7/2)))*(1-(3/2)*sin(i)^2));
a = fzero(a_fun,a0);
end