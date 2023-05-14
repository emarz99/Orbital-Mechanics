%% Assignemnt 2 - Planetary explorer mission

clear
close all
clc 

%% To call data
dataClean2;
kep0 = [a, e, i, OM, om, th0];

%% Nominal Orbit and nominal ground tracks

T = 2*pi*sqrt(a^3/mu_E);  % Orbital period [1/s]
theta_G0 = 0;             % True anomaly of the Greenwich meridian [rad]
t0 = 0;                   % Initial time [s]
k1 = 1;                   % One revolution of the satllite [-]
k2 = (3600*24)/T;         % Revolution of the satellite in 1 day [-]
k3 = 10*k2;               % Revolution of the satellite in 10 days [-]
N = 10000;                % Number of points considered [-]

% To propagate the nominal orbit and to compute and plot the nominal ground track
[r, v, lon, lat] = groundTrack(kep0, mu_E, t0, om_E, theta_G0, k1, N);

[~, ~, lon1, lat1] = groundTrack(kep0, mu_E, t0, om_E, theta_G0, k2, N);

[~, ~, lon2, lat2] = groundTrack(kep0, mu_E, t0, om_E, theta_G0, k3, N);

% To plot the nominal orbit
Terra_3D;
hold on
plot3(r(1,:),r(2,:),r(3,:));
title('Nominal Orbit')

%% Repeated Ground Track

% To compute the modified semi-major axis for the repeating ground track
a_rep = aRepGroundTrack(ratio,om_E,mu_E);
a_repJ2 = aRepGroundTrack_J2(ratio,om_E,mu_E,J2,R,e,i,a);
kep1 = kep0;
kep1(1) = a_rep;
k4 = 7;

% To propagate the orbit with repeating ground track and to compute and plot the repeating ground track
[r_rep, v_rep, lon_r, lat_r] = groundTrack(kep1, mu_E, t0, om_E, theta_G0, k4, N);

[~, ~, lon1_r, lat1_r] = groundTrack(kep1, mu_E, t0, om_E, theta_G0, k2, N);

[~, ~, lon2_r, lat2_r] = groundTrack(kep1, mu_E, t0, om_E, theta_G0, k3, N);



%% Perturbed ground tracks

% To propagate the orbit under perturbations and to compute and plot the ground track
[~, ~, lon_p, lat_p] = pertGroundTrack(kep0, mu_E, t0, om_E, theta_G0, k1, N, J2, R, A_M, cD);

[~, ~, lon1_p, lat1_p] = pertGroundTrack(kep0, mu_E, t0, om_E, theta_G0, k2, N, J2, R, A_M, cD);

[~, ~, lon2_p, lat2_p] = pertGroundTrack(kep0, mu_E, t0, om_E, theta_G0, k3, N, J2, R, A_M, cD);


% To propagate the orbit under perturbations and to compute and plot the repeating ground track
[~, ~, lon_rp, lat_rp] = pertGroundTrack(kep1, mu_E, t0, om_E, theta_G0, k4, N, J2, R, A_M, cD);

[~, ~, lon1_rp, lat1_rp] = pertGroundTrack(kep1, mu_E, t0, om_E, theta_G0, k2, N, J2, R, A_M, cD);

[~, ~, lon2_rp, lat2_rp] = pertGroundTrack(kep1, mu_E, t0, om_E, theta_G0, k3, N, J2, R, A_M, cD);

%% To plot the ground tracks

fig = imread('planisfero.jpg');

figure()
imagesc([-180 180], [90 -90], fig)
hold on
set(gca,'ydir','normal');
plot(lon,lat,'g');
plot(lon_p, lat_p, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon(1),lat(1),80,'r','filled');
scatter(lon(end),lat(end),80,'g','filled');
scatter(lon_p(end),lat_p(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-90 90])
grid on
title('Nominal and perturbed ground track - 1 orbit')
legend('Nominal Orbit Ground Track', 'Perturbed ground track','Initial point','Unperturbed final point', 'Perturbed final point')

% Zoom 
figure()
plot(lon,lat,'b');
hold on
plot(lon_p, lat_p, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon(1),lat(1),80,'r','filled');
scatter(lon(end),lat(end),80,'b','filled');
scatter(lon_p(end),lat_p(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-8 10])
grid on
legend('Nominal Orbit Ground Track', 'Perturbed ground track','Initial point','Unperturbed final point', 'Perturbed final point')

figure()
imagesc([-180 180], [90 -90], fig)
hold on
set(gca,'ydir','normal');
plot(lon1,lat1,'g');
plot(lon1_p, lat1_p, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon1(1),lat1(1),80,'r','filled');
scatter(lon1(end),lat1(end),80,'g','filled');
scatter(lon1_p(end),lat1_p(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-90 90])
grid on
title('Nominal and perturbed ground track - 1 day')
legend('Nominal Orbit Ground Track', 'Perturbed ground track','Initial point','Unperturbed final point', 'Perturbed final point')

% Zoom 
figure()
plot(lon1,lat1,'b');
hold on
plot(lon1_p, lat1_p, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon1(1),lat1(1),80,'r','filled');
scatter(lon1(end),lat1(end),80,'b','filled');
scatter(lon1_p(end),lat1_p(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-8 10])
grid on
legend('Nominal Orbit Ground Track', 'Perturbed ground track','Initial point','Unperturbed final point', 'Perturbed final point')

figure()
imagesc([-180 180], [90 -90], fig)
hold on
set(gca,'ydir','normal');
plot(lon2,lat2,'g');
plot(lon2_p, lat2_p, '--','Color',[0.9290 0.6940 0.1250])
scatter(lon2(1),lat2(1),80,'r','filled');
scatter(lon2(end),lat2(end),80,'g','filled');
scatter(lon2_p(end),lat2_p(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-90 90])
grid on
title('Nominal and perturbed ground track - 10 days')
legend('Nominal Orbit Ground Track', 'Perturbed ground track','Initial point','Unperturbed final point', 'Perturbed final point')

% Zoom 
figure()
plot(lon2,lat2,'b');
hold on
plot(lon2_p, lat2_p, 'Color',[0.9290 0.6940 0.1250])
scatter(lon2(1),lat2(1),80,'r','filled');
scatter(lon2(end),lat2(end),80,'b','filled');
scatter(lon2_p(end),lat2_p(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-8 10])
grid on
legend('Nominal Orbit Ground Track', 'Perturbed ground track','Initial point','Unperturbed final point', 'Perturbed final point')
%%
figure()
imagesc([-180 180], [90 -90], fig)
hold on
set(gca,'ydir','normal');
plot(lon_r,lat_r,'g');
plot(lon_rp, lat_rp, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon_r(1),lat_r(1),80,'r','filled');
scatter(lon_r(end),lat_r(end),80,'g','filled');
scatter(lon_rp(end),lat_rp(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-90 90])
grid on
title('Nominal and perturbed repeating ground track')
legend('Unperturbed Repeating Ground Track', 'Perturbed repeating ground track','Initial point','Unperturbed final point', 'Perturbed final point')

% Zoom 
figure()
plot(lon_r,lat_r,'b');
hold on
plot(lon_rp, lat_rp, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon_r(1),lat_r(1),80,'r','filled');
scatter(lon_r(end),lat_r(end),80,'b','filled');
scatter(lon_rp(end),lat_rp(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-8 10])
grid on
legend('Unperturbed Repeating Ground Track', 'Perturbed repeating ground track','Initial point','Unperturbed final point', 'Perturbed final point')

figure()
imagesc([-180 180], [90 -90], fig)
hold on
set(gca,'ydir','normal');
plot(lon1_r,lat1_r,'g');
plot(lon1_rp, lat1_rp, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon1_r(1),lat1_r(1),80,'r','filled');
scatter(lon1_r(end),lat1_r(end),80,'g','filled');
scatter(lon1_rp(end),lat1_rp(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-90 90])
grid on
title('Nominal and perturbed repeating ground track - 1 day')
legend('Unperturbed Repeating Ground Track', 'Perturbed repeating ground track','Initial point','Unperturbed final point', 'Perturbed final point')

% Zoom 
figure()
plot(lon1_r,lat1_r,'b');
hold on
plot(lon1_rp, lat1_rp, '--','Color', [0.9290 0.6940 0.1250])
scatter(lon1_r(1),lat1_r(1),80,'r','filled');
scatter(lon1_r(end),lat1_r(end),80,'b','filled');
scatter(lon1_rp(end),lat1_rp(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-8 10])
grid on
legend('Unperturbed Repeating Ground Track', 'Perturbed repeating ground track','Initial point','Unperturbed final point', 'Perturbed final point')

figure()
imagesc([-180 180], [90 -90], fig)
hold on
set(gca,'ydir','normal');
plot(lon2_r,lat2_r,'g');
plot(lon2_rp, lat2_rp,'Color', [0.9290 0.6940 0.1250])
scatter(lon2_r(1),lat2_r(1),80,'r','filled');
scatter(lon2_r(end),lat2_r(end),80,'g','filled');
scatter(lon2_rp(end),lat2_rp(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-90 90])
grid on
title('Nominal and perturbed repeating ground track - 10 days')
legend('Unperturbed Repeating Ground Track', 'Perturbed repeating ground track','Initial point','Unperturbed final point', 'Perturbed final point')

% Zoom 
figure()
plot(lon2_r,lat2_r,'b');
hold on
plot(lon2_rp, lat2_rp,'Color', [0.9290 0.6940 0.1250])
scatter(lon2_r(1),lat2_r(1),80,'r','filled');
scatter(lon2_r(end),lat2_r(end),80,'b','filled');
scatter(lon2_rp(end),lat2_rp(end),80,[0.9290 0.6940 0.1250],'filled');
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-8 10])
grid on
legend('Unperturbed Repeating Ground Track', 'Perturbed repeating ground track','Initial point','Unperturbed final point', 'Perturbed final point')


%% Perturbed Orbit Propagation
k = 50;                                                   % number of revolution of the satellite [-]
N = 10000;
tspan = linspace(t0,k*T,N);                               % [s]
options = odeset('RelTol',1e-13,'AbsTol',1e-14);

tic;
% To solve the ordinary differential equation in Keplerian elements with respect to time with Gauss planetary equation
[time, kep_pert] = odeSolver(kep0, mu_E, J2, R, A_M, cD, om_E, tspan, options);
toc

r_Gauss = zeros(3, N);
v_Gauss = zeros(3, N);
r_norm = zeros(1,N);
for mm = 1:length(time)
    [r_Gauss(:,mm), v_Gauss(:,mm)] = kep2car(kep_pert(mm,1), kep_pert(mm,2), kep_pert(mm,3), kep_pert(mm,4), kep_pert(mm,5), kep_pert(mm,6), mu_E);
    r_norm(mm) = norm(r_Gauss(:,mm));
end

tic;
% To solve the ordinary differential equation in Cartesian coordinates
[r0,v0] = kep2car(a,e,i,OM,om,th0,mu_E);
[~,state] = ode113(@(t,s) tbp_pert_ode(t, s, mu_E, a_Per_cart(t, s, mu_E, J2, R, A_M, cD, om_E)), tspan, [r0;v0], options);
toc
r_Cart = [state(:,1),state(:,2),state(:,3)];
v_Cart = [state(:,4),state(:,5),state(:,6)];


colourRGB = jet(k);                     %Generating colours to be used using jet colormap
colourtimes = N/k;                      %Determining num of times each colour will be used
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
colormap(colourRGB);                    %using the custom colormap of the colors we used
grid on; 
caxis([1  k]);

% To plot the solution in Cartesian coordinates
Terra_3D
hold on
scatter3(r_Cart(:,1),r_Cart(:,2),r_Cart(:,3), 2, colourRGB(colourind,:), 'filled')
title('Orbit propagation with Cartesian cooordinates')
colormap(colourRGB);                    %using the custom colormap of the colors we used
grid on; 
caxis([1  k]);


%% Evolution of Keplerian Elements
% Keplerian elements from Gauss method
a_G = kep_pert(:,1)';
e_G = kep_pert(:,2)';
i_G = kep_pert(:,3)';
OM_G = kep_pert(:,4)';
om_G = kep_pert(:,5)';
th_G = kep_pert(:,6)';

% Keplerian elements from Cartesian method
a_C = zeros(1, N);
e_C = zeros(1, N);
i_C = zeros(1, N);
OM_C = zeros(1, N);
om_C = zeros(1, N);
th_C = zeros(1, N);
for m = 1:length(time)
    [a_C(m), e_C(m), i_C(m), OM_C(m), om_C(m), th_C(m)] = car2kep(r_Cart(m,:)', v_Cart(m,:)', mu_E);
end
th_C = unwrap(th_C);


%% To evaluate the error between the two methods
T_plot=linspace(0, k, N);         % Periods

err_a = abs(a_C-a_G)/a;           % error on the semi-major axis
err_e = abs(e_C-e_G);             % error on the eccentricity
err_i = abs(i_C-i_G)/(2*pi);      % error on the inclination
err_OM = abs(OM_C-OM_G)/(2*pi);   % error on the RAAN
err_om = abs(om_C-om_G)/(2*pi);   % error on the argument of pericentre
err_th = abs((th_G-th_C)./th_G);  % error on the true anomaly

% To plot the errors
figure
subplot(2,3,1)
semilogy(T_plot,err_a)
title('Error on the semi-major axis')
xlabel('T -- Periods')
ylabel('\mida_C-a_G\mid/a_0')
grid on

hold on
subplot(2,3,2)
semilogy(T_plot,err_e)
title('Error on the eccentricity')
xlabel('T -- Periods')
ylabel('\mide_C-e_G\mid')
grid on

hold on
subplot(2,3,3)
semilogy(T_plot,err_i)
title('Error on the inclination')
xlabel('T -- Periods')
ylabel('\midi_C-i_G\mid/2\pi')
grid on

hold on
subplot(2,3,4)
semilogy(T_plot,err_OM)
title('Error on the RAAN')
xlabel('T -- Periods')
ylabel('\mid\Omega_C-\Omega_G\mid/2\pi')
grid on

hold on
subplot(2,3,5)
semilogy(T_plot,err_om)
title('Error on the argument of pericentre')
xlabel('T -- Periods')
ylabel('\mid\omega_C-\omega_G\mid/2\pi')
grid on

hold on
subplot(2,3,6)
semilogy(T_plot,err_th)
title('Error on the true anomaly')
xlabel('T -- Periods')
ylabel('\mid\theta_C-\theta_G\mid/\theta_G')
grid on


%% Filters
[f_co_a] = lowPassFilter(a_G);
a_lf = movmean(a_G,f_co_a);
a_lf1 = movmean(a_lf,2*f_co_a);

[f_co_e] = lowPassFilter(e_G);
e_lf = movmean(e_G,f_co_e);
e_lf1 = movmean(e_lf,2*f_co_e);

[f_co_i] = lowPassFilter(i_G);
i_lf = movmean(i_G,f_co_i);
i_lf1 = movmean(i_lf,2*f_co_i);

[f_co_OM,f_co_th] = lowPassFilterMono(kep_pert, mu_E, R, J2, A_M, cD, om_E);
OM_lf = movmean(OM_G,f_co_OM);
OM_lf1 = movmean(OM_lf,2*f_co_OM);

th_lf = movmean(th_G,f_co_th);
th_lf1 = movmean(th_lf,2*f_co_th);

[f_co_om] = lowPassFilter(om_G);
om_lf = movmean(om_G,f_co_om);
om_lf1 = movmean(om_G,2*f_co_om);

figure
plot(T_plot,a_G,T_plot,a_lf,T_plot,a_lf1)
title('Semi-major axis evolution')
legend('Short-period effect','Long-period effect','Secular effect')
ylabel('a [km]')
xlabel('T -- Periods')
grid on
axes('position',[.60 .35 .25 .25])
box on 
indexOfInterest = (T_plot < 30) & (T_plot > 26); 
plot(T_plot(indexOfInterest),a_G(indexOfInterest),T_plot(indexOfInterest),a_lf(indexOfInterest),T_plot(indexOfInterest),a_lf1(indexOfInterest))
axis tight

figure
plot(T_plot,e_G, T_plot, e_lf, T_plot, e_lf1)
title('Eccentricity evolution')
legend('Short-period effect','Long-period effect','Secular effect')
ylabel('e [-]')
xlabel('T -- Periods')
grid on
axes('position',[.60 .35 .25 .25])
box on 
indexOfInterest = (T_plot < 30) & (T_plot > 26); 
plot(T_plot(indexOfInterest),e_G(indexOfInterest),T_plot(indexOfInterest),e_lf(indexOfInterest),T_plot(indexOfInterest),e_lf1(indexOfInterest))
axis tight

figure
plot(T_plot,rad2deg(i_G), T_plot, rad2deg(i_lf), T_plot, rad2deg(i_lf1))
title('Inclination evolution')
legend('Short-period effect','Long-period effect','Secular effect')
ylabel('i [deg]')
xlabel('T -- Periods')
grid on
axes('position',[.60 .20 .25 .25])
box on 
indexOfInterest = (T_plot < 30) & (T_plot > 26); 
plot(T_plot(indexOfInterest),i_G(indexOfInterest),T_plot(indexOfInterest),i_lf(indexOfInterest),T_plot(indexOfInterest),i_lf1(indexOfInterest))
axis tight

figure
plot(T_plot, rad2deg(OM_G), T_plot, rad2deg(OM_lf), T_plot, rad2deg(OM_lf1))
title('Right ascension of the ascending node evolution')
legend('Short-period effect','Long-period effect','Secular effect')
ylabel('\Omega [deg]')
xlabel('T -- Periods')
grid on
axes('position',[.30 .20 .25 .25])
box on 
indexOfInterest = (T_plot < 30) & (T_plot > 26); 
plot(T_plot(indexOfInterest),rad2deg(OM_G(indexOfInterest)),T_plot(indexOfInterest),rad2deg(OM_lf(indexOfInterest)),T_plot(indexOfInterest),rad2deg(OM_lf1(indexOfInterest)))
axis tight

figure
plot(T_plot,rad2deg(om_G), T_plot, rad2deg(om_lf), T_plot, rad2deg(om_lf1))
title('Argument of pericentre evolution')
legend('Short-period effect','Long-period effect','Secular effect')
ylabel('\omega [deg]')
xlabel('T -- Periods')
grid on
axes('position',[.60 .20 .25 .25])
box on 
indexOfInterest = (T_plot < 30) & (T_plot > 26); 
plot(T_plot(indexOfInterest),rad2deg(om_G(indexOfInterest)),T_plot(indexOfInterest),rad2deg(om_lf(indexOfInterest)),T_plot(indexOfInterest),rad2deg(om_lf1(indexOfInterest)))
axis tight


figure
plot(T_plot,rad2deg(th_G), T_plot, rad2deg(th_lf), T_plot, rad2deg(th_lf1))
title('True anomaly evolution')
legend('Short-period effect','Long-period effect','Secular effect')
ylabel('\theta [deg]')
xlabel('T -- Periods')
grid on
axes('position',[.60 .20 .25 .25])
box on 
indexOfInterest = (T_plot < 30) & (T_plot > 26); 
plot(T_plot(indexOfInterest),rad2deg(th_G(indexOfInterest)),T_plot(indexOfInterest),rad2deg(th_lf(indexOfInterest)),T_plot(indexOfInterest),rad2deg(th_lf1(indexOfInterest)))
axis tight



%% Comparison of the model with real data
% DATA LOADED from 2017-9-4 to 2017-11-4 ----> circa 140 revolutions
load('Horizon_data.mat');

initialdate=date2mjd2000([2017 09 04 0 0 0]);
finaldate=date2mjd2000([2017 11 04 0 0 0]);

a_eph=str2double(Horizon_data(:,14));
e_eph=str2double(Horizon_data(:,5));
i_eph=str2double(Horizon_data(:,7));
OM_eph=str2double(Horizon_data(:,8));
om_eph=str2double(Horizon_data(:,9));
th_eph=str2double(Horizon_data(:,13));

kep_eph=[a_eph, e_eph, i_eph, OM_eph, om_eph, th_eph];

N_eph=length(kep_eph);
tend = (finaldate-initialdate)*24*3600;          % Duration of the time interval [s]
k_eph = tend/(2*pi*sqrt(a_eph(1)^3/mu_E));        % Number of orbits of the satellite

tspan = linspace(t0, tend, N_eph);               % Time interval with N points (for the model) [s]
T_plot2=linspace(0, k_eph, N_eph);



%% Perturbed Orbit Propagation
kep0 = kep_eph(1,:);
kep0(3:6) = deg2rad(kep0(3:6));
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Orbit propagation with Gauss Planetary equations
[time, kep_pert1] = odeSolver(kep0, mu_E, J2, R, A_M, cD, om_E, tspan, options);

% Orbit propagation in Cartesian coordinates
[r0,v0] = kep2car(kep0(1),kep0(2),kep0(3),kep0(4),kep0(5),kep0(6),mu_E);
[~,state1] = ode113(@(t,s) tbp_pert_ode(t, s, mu_E, a_Per_cart(t, s, mu_E, J2, R, A_M, cD, om_E)), tspan, [r0;v0], options);
r_Cart1 = [state1(:,1),state1(:,2),state1(:,3)];
v_Cart1 = [state1(:,4),state1(:,5),state1(:,6)];


%% Evolution of Keplerian Elements

% Keplerian elements from Gauss method
a_G1 = kep_pert1(:,1);
e_G1 = kep_pert1(:,2);
i_G1 = kep_pert1(:,3);
OM_G1 = kep_pert1(:,4);
om_G1 = kep_pert1(:,5);
th_G1 = wrapTo2Pi(kep_pert1(:,6));


% Keplerian elements from Cartesian method
a_C1 = zeros(1, N_eph);
e_C1 = zeros(1, N_eph);
i_C1 = zeros(1, N_eph);
OM_C1 = zeros(1, N_eph);
om_C1 = zeros(1, N_eph);
th_C1 = zeros(1, N_eph);
for m = 1:length(tspan)
    [a_C1(m), e_C1(m), i_C1(m), OM_C1(m), om_C1(m), th_C1(m)] = car2kep(r_Cart1(m,:)', v_Cart1(m,:)', mu_E);
end



% Comparison of the results 
figure
plot(T_plot2, a_G1, T_plot2, a_C1,'--', T_plot2, kep_eph(:,1))
title('Semi-major axis')
legend('Gauss Equations','Cartesian Coordinates', 'Real evolution')
ylabel('a [km]')
xlabel('T -- Periods')
grid on

figure
plot(T_plot2, e_G1, T_plot2, e_C1,'--', T_plot2, kep_eph(:,2))
title('Eccentricity')
legend('Gauss Equations','Cartesian Coordinates', 'Real evolution')
ylabel('e [-]')
xlabel('T -- Periods')
grid on

figure
plot(T_plot2, rad2deg(i_G1), T_plot2, rad2deg(i_C1),'--', T_plot2, kep_eph(:,3))
title('Inclination')
legend('Gauss Equations','Cartesian Coordinates', 'Real evolution')
ylabel('i [deg]')
xlabel('T -- Periods')
grid on

figure
plot(T_plot2, rad2deg(OM_G1), T_plot2, rad2deg(OM_C1),'--', T_plot2, kep_eph(:,4))
title('Right Ascension of Ascending Node')
legend('Gauss Equations','Cartesian Coordinates', 'Real evolution')
ylabel('\Omega [deg]')
xlabel('T -- Periods')
grid on

figure
plot(T_plot2, rad2deg(om_G1), T_plot2, rad2deg(om_C1),'--', T_plot2, kep_eph(:,5))
title('Argument of the pericentre')
legend('Gauss Equations','Cartesian Coordinates', 'Real evolution')
ylabel('\omega [deg]')
xlabel('T -- Periods')
grid on

figure
plot(T_plot2, rad2deg(th_G1), T_plot2, rad2deg(th_C1),'--', T_plot2, kep_eph(:,6))
title('True anomaly')
legend('Gauss Equations','Cartesian Coordinates', 'Real evolution')
ylabel('\theta [deg]')
xlabel('T -- Periods')
grid on