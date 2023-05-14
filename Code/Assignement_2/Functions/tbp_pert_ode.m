function ds = tbp_pert_ode(~, s, mu, ap_cart)

% tbp_pert_ode - Differential equation for the perturbed two-body problem
% in Cartesian reference frame
% Considered perturbatiosn: zonal harmonics (two dimensions) and atmospheric drag

% INPUT:
% s       [1x6] State vector in cartesian coordinates
% mu      [1x1] Gravitational parameter [km^3/s^2]
% ap_cart [3x1] Perturbing acceleration in cartesian coordinates [km/s^2]

% OUTPUT:
% ds      [1x6] Derivative of the state vector 

% AUTHORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso


% To extract the position and the velocity vectors from the state vector in
% Cartesian reference frame

r = s(1:3); % [km]
v = s(4:6); % [km/s]

r_norm = norm(r); % Absolute value of the position vector [km] 

% To compute the derivative of the state vector in Cartesian reference
% frame

ds = [v(1); v(2); v(3); -mu/r_norm^3*r(1) + ap_cart(1); -mu/r_norm^3*r(2) + ap_cart(2); -mu/r_norm^3*r(3) + ap_cart(3)];

end
