function [r, v, lon, lat] = groundTrack(kep0, mu, t0, om_E, theta_G0, k, N)

% groundTrack: propagate the orbit and compute latitude and longitude 

% INPUT:
% kep0       [1x6] Vector of initial keplerian elements
% mu         [1x1] Gravitational parameter [km^3/s^2]
% t0         [1x1] Initial time [s]
% om_E       [1x1] Earth's rotation velocity [rad/s]
% theta_G0   [1x1] 
% k          [1x1] Number of satellite orbits
% N          [1x1] Number of points considered for the integration
% 
% OUTPUT:
% r          [3xN] Position vector [km]
% v          [3xN] Velocity vector [km/s]
% lon        [Nx1] Longitude [deg]
% lat        [Nx1] Latitude [deg]

% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso


% Orbit propagation

a = kep0(1);   % Semi-major axis [km]
e = kep0(2);   % Eccentricity [-]
i = kep0(3);   % Inclination [rad]
OM = kep0(4);  % RAAN [rad]
om = kep0(5);  % Argument of pericenter [rad]
th0 = kep0(6); % True Anomaly [rad]

T = 2*pi*sqrt(a^3/mu);      % Orbital period [1/s]

[r0,v0] = kep2car(a,e,i,OM,om,th0,mu);
tspan = linspace(0,k*T,N);

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[time,state] = ode45(@(t,s) tbp_ode(t,s,mu),tspan,[r0;v0], options);

r = [state(:,1),state(:,2),state(:,3)]'; 
v = [state(:,4),state(:,5),state(:,6)]'; 

% Conversion to RA and declination

r_norm = sqrt(dot(r,r));

delta = zeros(size(time,1),1);
alpha = zeros(size(time,1),1);
for ii = 1:size(r,2)
    delta(ii) = asin(r(3,ii)./r_norm(ii)); % [rad]
    if r(2,ii)/r_norm(ii) >0 
        alpha(ii) = acos((r(1,ii)./r_norm(ii))./cos(delta(ii))); % [rad]
    else
        alpha(ii) = 2*pi-acos((r(1,ii)./r_norm(ii))./cos(delta(ii)));
    end
end

% Conversion to longitude and latitude

theta_G = theta_G0+om_E*(time-t0); % [deg]
lon = wrapTo180(rad2deg(alpha)-rad2deg(theta_G)); % [deg]
lat = rad2deg(delta); % [deg]

% To delete lines

for g = 1:length(lon)-1
    gg = g+1;
    L = abs(lon(gg)-lon(g));
    if L > 350
        lon(gg) = NaN;
    end
end

end
