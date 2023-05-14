function [r, v, lon, lat] = pertGroundTrack(kep0, mu, t0, om_E, theta_G0, k, N, J2, R, A_M, cD)

% pertGroundTrack: propagate the perturbed orbit and compute latitude and
% longitude

% INPUT:
% kep0       [1x6] Vector of initial keplerian elements
% mu         [1x1] Gravitational parameter [km^3/s^2]
% t0         [1x1] Initial time [s]
% om_E       [1x1] Earth's rotation velocity [rad/s]
% theta_G0   [1x1] True anomaly od Greenwich meridian at time t0 [rad]
% k          [1x1] Number of satellite orbits
% N          [1x1] Number of points considered for the integration
% J2         [1x1] Zonal harmonics (two dimensions) [-]
% R          [1x1] Radius of the planet [km]
% A_M        [1x1] Area to mass ratio [km^2/kg]
% cD         [1x1] Drag coefficient [-]
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
% Guidetto Tommaso


% Orbit propagation

a = kep0(1);                % Semi-major axis [km]
T = 2*pi*sqrt(a^3/mu);      % Orbital period [1/s]
tspan = linspace(0,k*T,N);

options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[time, kep] = odeSolver(kep0, mu, J2, R, A_M, cD, om_E, tspan, options);

for mm = 1:length(time)
    [r(:,mm), v(:,mm)] = kep2car(kep(mm,1), kep(mm,2), kep(mm,3), kep(mm,4), kep(mm,5), kep(mm,6), mu);
end

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

theta_G = theta_G0+om_E*(time-t0);                 % True anomaly of Greenwich Meridian [rad]
lon = wrapTo180(rad2deg(alpha)-rad2deg(theta_G));  % Longitude [deg]
lat = rad2deg(delta);                              % Latitude [deg]

% To delete the lines

for g = 1:length(lon)-1
    gg = g+1;
    L = abs(lon(gg)-lon(g));
    if L > 350
        lon(gg) = NaN;
    end
end

end