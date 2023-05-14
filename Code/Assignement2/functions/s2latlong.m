function [alpha, delta, lon, lat] =  s2latlong(s, t, Gr_lon, w_planet)
%--------------------------------------------------------------------------
%   Calculates latitude, longitude, right ascension and declination from
%   the cartesian state of the satellite.
%--------------------------------------------------------------------------
%   Form:
%   [alpha, delta, lon, lat] =  s2latlong(s, t, Gr_lon, w_planet)
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   s             (:,6)   Cartesian state of the satellite (x, y, z, vx,
%                         vy, vz) [L, L/T][Km, Km/s]
%   t             (:,1)   Time vector corresponding to s. [T][s]
%   Gr_lon        [1]     Longitude of Greenwich meridian at t(1). [][rad]
%   w_planet      [1]     Angular velocity of the central planet[1/T][rad/s].
%  
%   -------
%   Outputs
%   -------
%   alpha         (:,1)   
%
%--------------------------------------------------------------------------
% Programmed by: Jaime Fernández Diz
%
% Date:                  20/12/2022
% Revision:              
% Tested by:
%--------------------------------------------------------------------------
alpha = zeros(1,size(s,1));
delta = zeros(1,size(s,1));
lon = zeros(1, size(s,1));

for i = 1:size(s,1)
    alpha(i) = atan2(s(i,2), s(i,1));
    delta(i) = asin(s(i,3)/norm(s(i,1:3)));
    lon(i) = mod(pi + alpha(i) - (Gr_lon + w_planet*t(i)), 2*pi) - pi;
end

lat = delta;

end