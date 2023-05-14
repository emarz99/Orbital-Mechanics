function a = a_J2(t, s)

% Provides the acceleration due to the J2 term of the geopotential
%
% PROTOTYPE:
% dy = a_J2(t, cartesian_state)
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% s[6x1] Cartesian state of the body (rx , ry , rz , vx , vy , vz ) [ L,
% L/T] [Km, Km/s]
% OUTPUT:
% dy [3x1] Acceleration due to J2 [L/T^2] [Km/s^2]
%
% CONTRIBUTORS:
% Jaime Fernández Diz
%
% VERSIONS
% 2022 12 22: First version

r = norm(s(1:3));
mu_earth = astroConstants(13);
R_earth = astroConstants(23);
J2 = astroConstants(9);

a = 3*J2*mu_earth*R_earth^2/(2*r^4)*[s(1)*(5*s(3)^2/r^2-1)/r, ...
                                     s(2)*(5*s(3)^2/r^2-1)/r, ...
                                     s(3)*(5*s(3)^2/r^2-3)/r]';

end