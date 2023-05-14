function [sdot] = TwoBodyPerturbed(t, s, mu, perturbations)
%--------------------------------------------------------------------------
%   Newton equation for two body problem with a given perturbation
%   acceleration. The perturbation should be a function of s and t.
%--------------------------------------------------------------------------
%   Form:
%   sdot = TwoBodyPerturbed(t_perturbed, s_perturbed, mu_earth, @(t, s) a_J2(t, s) + a_moon(t,s,date0)
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   t             [1]     Instant of time of the propagation. [T][s]
%   s             (1,6)   Cartesian state (x, y, z, vx, vy, vz). 
%                         [L, L, L, L/T, L/T, L/T][Km,Km,Km,Km/s,Km/s,Km/s].
%   mu            [1]     Gravitational parameter of the central body. 
%                         [L^3/M T^2][km^3/(kg*s^2)].
%   perturbations [1]     Pointer to a function which returns the
%                         perturbation acceleration in cartesian 
%                         coordinates.         
%
%   -------
%   Outputs
%   -------
%   sdot         (1,6)   Derivative of state vector(vx, vy, vz, ax, ay, az). 
%                        [L/T, L/T, L/T, L/T^2, L/T^2, L/T^2]
%                        [Km/s,Km/s,Km/s,Km/s^2,Km/s^2,Km/s^2].
%
%--------------------------------------------------------------------------
% Programmed by: Jaime Fernández Diz
%
% Date:                  22/12/2022
% Revision:              
% Tested by:
%--------------------------------------------------------------------------
r = norm(s(1:3));
sdot = [s(4),...
        s(5),...
        s(6),...
        -mu*s(1)/r^3,...
        -mu*s(2)/r^3,...
        -mu*s(3)/r^3]';

sdot(4:6) =  sdot(4:6) + perturbations(t, s);

end