function ap_rsw = a_Per_rsw(~, kep, mu, R, J2, A_M, cD, om_E)

%a_Per_rsw: compute the perturbing acceleration in radial-tangential-out of
%       plane reference frame. Perturbations considered are zonal harmonics
%       and atmospheric drag.

% INPUT:
% kep        [1x6] Vector of keplerian elements
% mu         [1x1] Gravitational parameter of the planet [km^3/s^2]
% R          [1x1] Radius of the planet [km]
% J2         [1x1] Zonal harmonics (two dimensions) [-]
% A_M        [1x1] Area to mass ratio [km^2/kg]
% cD         [1x1] Drag coefficient [-]
% om_E       [1x1] Earth's rotation velocity [rad/s]

% OUTPUT:
% ap         [1x3] Perturbing acceleration in rsw frame [km/s^2]

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetto Tommaso
% Pasquariello Chiara

% From Keplerian Elements to Cartesian Coortinates

a = kep(1);   % Semi-major axis [km]
e = kep(2);   % Eccentricity [-]
i = kep(3);   % Inclination [rad]
OM = kep(4);  % RAAN [rad]
om = kep(5);  % Argument of pericenter [rad]
th = wrapTo2Pi(kep(6));  % True Anomaly [rad]

[r_vec,v_vec] = kep2car(a, e, i, OM, om, th, mu);

r = norm(r_vec); % Norm of position vector [km]
v = norm(v_vec); % Norm of velocity vector [km/s]

% Zonal harmonic perturbing acceleration in radial-transversal-out of plane reference
% frame

aJ2_rsw = (-3/2)*(J2*mu*R^2)/(r^4)*[1-3*sin(i)^2*sin(th+om)^2; sin(i)^2*sin(2*(th+om)); sin(2*i)*sin(th+om)];

% Atmospheric drag perturbing acceleration

aDrag_cart = drag(A_M, cD, om_E, r, R, r_vec, v_vec); % in cartesian frame
aDrag_rsw = cartesian2rsw(kep, aDrag_cart); % rotation to rsw frame

% Total perturbing acceleration

ap_rsw = aJ2_rsw + aDrag_rsw;

end

