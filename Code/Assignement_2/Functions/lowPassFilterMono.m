function [f_co_OM, f_co_th] = lowPassFilterMono(kep, mu, R, J2, A_M, cD, om_E)

% lowPassFilter - Select a proper cut-off frequency for a quasi-monothonic
% oscillations. In particular are evaluated for the RAAN and the true anomaly
% for the specific given orbit of the assignment

% INPUT:
% kep    [1x6] Keplerian elements vector
% mu     [1x1] Planetary constants of the planet [km^3/s^2]
% R      [1x1] Mean radius of the planet [km]
% J2     [1x1] Zonal harmonics (two dimensions) [-]
% A_M    [1x1] Cross area to mass ratio (for the drag computation) [m^2/kg]
% cD     [1x1] Drag coefficient (for the drag computation) [-]
% om_E   [1x1] Angular velocity of the planet along h [rad/s]

% OUTPUT:
% f_co_OM [1x1] cut-off frequency for the RAAN [Hz]
% f_co_th [1x1] cut-off frequency for the true anomaly [Hz]

% AUTHORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

% To evaluate the derivatives of the keplerian elements
KEP_dot = [];

for ii = 1:length(kep(:,1))
    a_p = a_Per_rsw(0, kep(ii,:), mu, R, J2, A_M, cD, om_E);
    [kep_dot] = gaussEoM_rsw( NaN, kep(ii,:), mu, a_p);
    KEP_dot = [KEP_dot,kep_dot];
end

% Filter for the RAAN

OM_dot = KEP_dot(4,:);

% Looking at OM_dot and OM evolutions in time it is possible to select a
% repeating value of OM_dot that is used to count how many points are in
% the intervall

J = [];
for kk = 1:length(kep(:,1))
    if OM_dot(kk) > -1e-9
        j = kk; % index of the element at which OM_dot cross the selected value
        J = [J;j];
    end
end

% To compute the number of points in each intervall

M = [];
for l = 1:(length(J)-1)
    m = J(l+1) - J(l);
    if m == 1
        m = NaN;
    end
    M = [M;m];
end

% To compute the cut-off frequency for the RAAN

f_co_OM = round(mean(M,'omitnan'));


% Filter for the true anomaly

th_dot = KEP_dot(6,:);

% Looking at th_dot and th evolutions in time it is possible to select a
% repeating value of th_dot that is used to count how many points are in
% the intervall

G = [];
for ff = 1:length(kep(:,1))
    if th_dot(ff) < 8e-5 % index of the element at which th_dot cross the selected value
        g = ff;
        G = [G;g];
    end
end

% To compute the number of points in each intervall

D = [];
for n = 1:(length(G)-1)
    d = G(n+1) - G(n);
    if d == 1
        d = NaN;
    end
    D = [D;d];
end

% To compute the cut-off frequency for the true anomaly

f_co_th = 2*round(mean(D,'omitnan'));
