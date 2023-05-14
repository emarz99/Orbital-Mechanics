function [kep_dot] = gaussEoM_rsw( ~, kep, mu, a_p)

% gaussEoM_tnh - Compute the Gauss planetary equations in the RSW reference 
%                frame starting from the Keplerian elements 

% INPUT:
% kep    [1x6] Keplerian elements:
% mu     [1x1] Planetary constant of the planet [km^3/s^2]
% a_p    [1x3] is the sum of J2 and drag accelertation vector

% OUTPUT:
% kep_dot [1x6] derivatives of the Keplerian elements

% AUTHORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

a = kep(1);  % Semi-major axis [km]
e = kep(2);  % Eccentricity [-]
i = kep(3);  % Inclination of the orbital plane [rad]
OM = kep(4); % Right ascension of the ascending node [rad]
om = kep(5); % Argument of pericentre [rad]
th = wrapTo2Pi(kep(6)); % True anomaly [rad]

% Compute the perturbing acceleration in the RSW frame

ar = a_p(1);
as = a_p(2);
aw = a_p(3);

p = a*(1-e^2);       % Semilatus rectum [km]
r = p/(1+e*cos(th)); % Absolute value of the position [km]
h = sqrt(p*mu);      % Absolute value of the angular momentum [km^2/s]

% Gauss planetary equations formula:

a_dot = (2*(a^2)/h)*(e*sin(th)*ar+p*as/r);
e_dot = (1/h)*(p*sin(th)*ar+((p+r)*cos(th)+r*e)*as);
i_dot = r*cos(th+om)*aw/h;
OM_dot = r*sin(th+om)*aw/(h*sin(i));
om_dot = ((1/(h*e))*(-p*cos(th)*ar+(p+r)*sin(th)*as))-(r*sin(th+om)*cos(i)*aw/(h*sin(i)));
th_dot = (h/(r^2))+(1/(e*h))*(p*cos(th)*ar-(p+r)*sin(th)*as);

% To save each darivative in a vector

kep_dot = [a_dot, e_dot, i_dot, OM_dot, om_dot, th_dot]';

end
