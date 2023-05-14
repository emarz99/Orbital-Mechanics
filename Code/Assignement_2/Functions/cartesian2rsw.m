function [v_rsw] = cartesian2rsw(kep, v_car)

% aRepGroundTrack_J2: allow to pass from a vector in cartesian-frame to the
%                      same vector expressed in rsw-frame

% INPUT:
% kep        [1x6] Vector of keplerian elements
% v_car       [3]  Vector in Cartesian frame

%OUTPUT:    
% V_rsw       [3]  Vector in RSW frame

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso
% Pasquariello Chiara

a = kep(1);     % Semi-major axis [km]
e = kep(2);     % Eccentricity [-]
i = kep(3);     % Inclination [rad]
OM = kep(4);    % RAAN [rad]
om = kep(5);    % Argument of pericenter [rad]
th = kep(6);    % True Anomaly [rad]

% Definition of the rotation matrices

R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1]; % Rotation of OM around third axis
R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; % Rotation of i around first axis
R3_omth = [cos(om+th) sin(om+th) 0; -sin(om+th) cos(om+th) 0; 0 0 1]; % Rotation of om around third axis

% Rotation of vector from rsw to ECEI frame
v_rsw = R3_omth*R1_i*R3_OM*v_car;

end