function ds = tbp_ode(~, s, mu)

% tbp_ode ODE system for the two-body problem

% INPUT:
% s       [1x6] State vector in cartesian coordinates
% mu      [1x1] Gravitational parameter [km^3/s^2]

% OUTPUT:
% ds      [1x6] Derivative of the state vector 

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso
% Pasquariello Chiara

r = s(1:3);       % Position vector [km]
v = s(4:6);       % Velocity vector [km/s]
r_norm = norm(r);

ds = [v(1); v(2); v(3); -mu/r_norm^3*r(1); -mu/r_norm^3*r(2); -mu/r_norm^3*r(3)];

end