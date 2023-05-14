function a = aRepGroundTrack(ratio, om_E, mu)

% aRepGroundTrack: compute the semi-major axis for repeating ground track

% INPUT:
% ratio   [1x1] Ratio between number of revolution of the satellite and 
%               number of rotations of the planet [-]
% om_E    [1x1] Earth's rotation velocity [rad/s]
% mu      [1x1] Gravitational parameter [km^3/s^2]

% OUTPUT:
% a       [1x1] Semi-major axis of modified orbit [km]

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso
% Pasquariello Chiara

n = om_E*(ratio);    % Mean motion of modified orbit [1/s]
a = (mu/n^2)^(1/3);  % Modified semi-major axis [km]

end