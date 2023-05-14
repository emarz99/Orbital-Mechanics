function ap = drag(A_M, cD, om_E, r, R, r_vec, v_vec)

% drag: compute the value of atmospheric drag acceleration along the
% direction of the relative velocity between satellite and air particles
% Cartesian coordinates

% INPUT:
% A_M        [1x1] Area to mass ratio [km^2/kg]
% cD         [1x1] Drag coefficient [-]
% om_E       [1x1] Earth's rotation velocity [rad/s]
% r          [1x1] Norm of th position vector [km]
% R          [1x1] Radius of the planet [km]
% r_vec      [3x1] Position vector [km]
% v_vec      [3x1] Velocity vector [km/s]

% OUTPUT:
% ap         [1x1] Atmospheric drag acceleration along t [km/s^2]

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso
% Pasquariello Chiara

% Computation of air density


alt = r - R;                  % Altitude [km]
if alt < 0
    disp('crash with earth')
    rho = nan;
else
    rho = expAtmModel(alt)*1e9;   % Density [kg/km^3]
end

% Vector of planet's rotational velocity

om_E_vec = [0; 0; om_E];      % [rad/s]

% Relative velocity between satellite and air particles

v_rel = (v_vec)-cross(om_E_vec,r_vec);  % Vector of relative velocity in Cartesian coordinates [km/s]
v_rel_norm = norm(v_rel); % Norm of relative velocity [km/s]

% Computation of air drag in Cartesian coordinates

ap = -0.5*A_M*cD*rho*v_rel_norm^2*(v_rel/v_rel_norm);

end
