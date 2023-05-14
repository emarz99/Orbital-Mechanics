function [ys,v] = kep2car(a,e,i,Omega,w,theta_0,mu)
% This function calculate sthe cartesian state vectors of position and
% velocity from the keplerian elements around a certain Celestial object
% with planetary constant 'mu'.
%
% PROTOTYPE:
%   [ys,v] = kep2car(a,e,i,Omega,w,theta_0,mu)
%
% INPUT:
%   a[1]        semi-major axis                                 [km]
%   e[1]        eccentricity                                    [-]
%   i[1]        inclination                                     [deg]
%   Omega[1]    right ascension of the ascending node (RAAN)    [deg]
%   w[1]        argument of pericenter                          [deg]
%   theta_0[1]  true anomaly                                    [deg]
%
% OUTPUT:
%   ys[3]       Position vector in Cartesian coordinates        [km]
%   v[3]        Velocity vector in Cartesian coordinates        [km/s]
%
% CONTRIBUTORS:
%     Carlos Albiñana
%
% VERSIONS:
%     2022-11-20:    First Version


if nargin < 6 % if less than 6 inputs, break
    error('At least 6 input arguments required.');
else
    if e >= 1 % if parabollic or hyperbolic orbit, break
        error('The eccentricity must be smaller than 1.');
    end
    

    % Convert angles in deg to rad
    i = i*pi/180;
    Omega = Omega*pi/180;
    w = w*pi/180;
    theta_0 = theta_0*pi/180;
    
    % Compute needed parameters for the transformation
    r = a*(1 - e^2)/(1 + e*cos(theta_0));
    rf = [r*cos(theta_0); r*sin(theta_0)];
    l1 = cos(Omega)*cos(w) - sin(Omega)*sin(w)*cos(i);
    l2 = -cos(Omega)*sin(w) - sin(Omega)*cos(w)*cos(i);
    m1 = sin(Omega)*cos(w) + cos(Omega)*sin(w)*cos(i);
    m2 = -sin(Omega)*sin(w) + cos(Omega)*cos(w)*cos(i);
    n1 = sin(w)*sin(i);
    n2 = cos(w)*sin(i);
    
    % Determine position vector
    X = [l1 l2; m1 m2; n1 n2]*rf;
    x = X(1);
    y = X(2);
    z = X(3);
    
    % Determine velocity vector
    H = sqrt(mu*a*(1 - e^2));
    v_x = mu/H*(-l1*sin(theta_0) + l2*(e + cos(theta_0)));
    v_y = mu/H*(-m1*sin(theta_0) + m2*(e + cos(theta_0)));
    v_z = mu/H*(-n1*sin(theta_0) + n2*(e + cos(theta_0)));    
end
v=[v_x,v_y,v_z];
ys=[x,y,z];