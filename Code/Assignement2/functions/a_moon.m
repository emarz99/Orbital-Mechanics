function a = a_moon(t, s, date)

% Provides the perturbation acceleration due to the Moon
%
% PROTOTYPE:
% dy = a_moon(t, cartesian_state, date)
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% s[6x1] Cartesian state of the body (rx , ry , rz , vx , vy , vz ) [ L,
% L/T] [Km, Km/s]
% date[1] MJD2000 of the beginning of the propagation [T] [Days]
% OUTPUT:
% dy [3x1] Acceleration due to Moon [L/T^2] [Km/s^2]
%
% CONTRIBUTORS:
% Jaime Fernández Diz
%
% VERSIONS
% 2022 12 22: First version

date = date + t/(24*3600);                    %T in seconds, date in MJD.
mu_moon = astroConstants(20);
[r_moon, ~] = ephMoon(date);
[axm, aym, azm] = thirdbody(mu_moon, s(1), s(2), s(3), r_moon(1), r_moon(2), r_moon(3));
a = [axm, aym, azm]';

end