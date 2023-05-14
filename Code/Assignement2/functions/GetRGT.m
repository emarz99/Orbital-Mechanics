function a = GetRGT(w_planet, mu_planet, k, m)
%--------------------------------------------------------------------------
%   Calculates the semimajor axis of a repeating ground track with given 
%   parameters.
%--------------------------------------------------------------------------
%   Form:
%   a = GetRGT(w_planet, mu_planet, k, m)
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   w_planet      [1]     Angular speed of the planet. [1/T][rad/s]
%   mu_planet     [1]     Gravitational parameter of the central body. 
%                         [L^3/M T^2][km^3/(kg*s^2)].
%   k             [1]     Number of complete revolutions of the satellite
%                         [#]
%   m             [1]     Number of complete revolutions of the planet
%   -------
%   Outputs
%   -------
%   a             [1]     SMA of the desired RGT. [L][Km]
%
%--------------------------------------------------------------------------
% Programmed by: Jaime Fernández Diz
%
% Date:                  20/12/2022
% Revision:              
% Tested by:
%--------------------------------------------------------------------------
T = m*2*pi/(k*w_planet);
a = (mu_planet*(T/(2*pi))^2)^(1/3);

end