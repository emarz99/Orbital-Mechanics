function ap_cart = a_Per_cart(~,state, mu, J2, R, A_M, cD, om_E)

%a_Per_cart: compute the perturbing acceleration in cartesian coordinates. 
%           Perturbations considered are zonal harmonics and atmospheric drag.

% INPUT:
% state      [1x6] State of the satellite
% mu         [1x1] Gravitational parameter of the planet [km^3/s^2]
% R          [1x1] Radius of the planet [km]
% J2         [1x1] Zonal harmonics (two dimensions) [-]
% A_M        [1x1] Area to mass ratio [km^2/kg]
% cD         [1x1] Drag coefficient [-]
% om_E       [1x1] Earth's rotation velocity [rad/s]

% OUTPUT:
% ap         [1x3] Perturbing acceleration in cartesian coordinates [km/s^2]

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetto Tommaso
% Pasquariello Chiara

r_vec = state(1:3);
v_vec = state(4:6);

r = norm(r_vec);

% Acceleration due to zonal harmonic
aJ2_cart = [((3*J2*mu*R^2*r_vec(1))/(2*r^5))*(5*(r_vec(3)^2)/r^2 -1);((3*J2*mu*R^2*r_vec(2))/(2*r^5))*(5*(r_vec(3)^2)/r^2 -1);((3*J2*mu*R^2*r_vec(3))/(2*r^5))*(5*(r_vec(3)^2)/r^2 -3)];

% Acceleration due to atmospheric drag in cartesian frame
aDrag_cart = drag(A_M, cD, om_E, r, R, r_vec, v_vec);

% Total perturbing acceleration
ap_cart = aJ2_cart + aDrag_cart;

end