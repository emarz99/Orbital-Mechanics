clear;
close all;
clc

%% To call Data
dataClean2;
kep0 = [a, e, i, OM, om, th0];
T = 2*pi*sqrt(a^3/mu_E);   % Orbital period [1/s]
theta_G0 = 0;              % True anomaly of the Greenwich meridian [rad]
t0 = 0;                    % Initial time [s]

J2=0;                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perturbed Orbit Propagation
k = 5000;                    % number of revolution of the satellite [-]
N = 500000;

% From Keplerian elements to Cartesian coordinates for the solution of the Gauss equation
[r_Gauss, v_Gauss, ~, ~, ~,] = pertGroundTrack(kep0, mu_E, t0, om_E, theta_G0, k, N, J2, R, A_M, cD);
close;

colourRGB = jet(k);                      %Generating colours to be used using jet colormap
colourtimes = N/k;                       %Determining num of times each colour will be used
colourind=[];
for ii=1:k
    cc=ii*ones(colourtimes, 1);
    colourind = [colourind, cc'];        %Linear indices of colours for newA
end

% To plot the solution with Gauss planetary equations
Terra_3D;
hold on
scatter3(r_Gauss(1,:),r_Gauss(2,:),r_Gauss(3,:), 2, colourRGB(colourind,:), 'filled')
title('Orbit propagation with Gauss planetary equations')
colormap(colourRGB);                      %using the custom colormap of the colors we used
grid on; 
caxis([1  k]);


%%  Verify Solution
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[time, kep_pert] = odeSolver(kep0, mu_E, J2, R, A_M, cD, om_E, tspan, options);
rp = kep_pert(:,1).*(1-kep_pert(:,2));
ra = kep_pert(:,1).*(1+kep_pert(:,2));
TK=linspace(0, k, N);

figure
plot(time, ra)
xlabel('T -- periods')
ylabel('r_a [Km]')
title('r_a VS revolutions')

figure
plot(time, rp)
xlabel('T -- periods')
ylabel('r_p [Km]')
title('r_p VS revolutions')