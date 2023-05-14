function sol = rp_solver(rp,v_inf_minus,v_inf_plus,delta_angle,muP)
% This function is used in combination with fzero to find the
% periapsis value of a powered fly-by given the incoming and outgoing
% planetocentric velocities and the angle between them.

% PROTOTYPE
%   sol = rp_solver(rp,v_inf_minus,v_inf_plus,delta_angle,mu)
%
% INPUT
%     rp[1]               Perigee - used as variable in fzero/fsolve [L]
%     v_inf_minus[1]      Incoming Planetocentric velocity magnitude [L/T]
%     v_inf_minus[1]      Incoming Planetocentric velocity magnitude [L/T]
%     delta_angle[1]      Turn angle formed between both planetocentric velocities [rad]
%     muP[1]              Gravitational Parameter of the primary [L^3/T^2]  
%     
% OUTPUT
%     sol[1]              Output of the equation that must become zero [rad]
%
% CONTRIBUTORS
%   Carlos Albiñana
%
% VERSIONS
% 08/12/2022: First Version

delta_minus = 2*asin(muP/(muP+rp*v_inf_minus^2));   % Turn Angle in incomming hyperbolic trajectory [rad]
delta_plus  = 2*asin(muP/(muP+rp*v_inf_plus^2));    % Turn Angle in outgoing hyperbolic trajectory  [rad]     

sol=0.5*(delta_minus+delta_plus)-delta_angle;       % Output of the function that must become zero  [rad]

end