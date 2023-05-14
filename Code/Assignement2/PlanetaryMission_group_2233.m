%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Assignment 2 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Planetary Explorer Mission

clear variables
close all
clc

addpath 'functions'  
addpath '../MATLAB functions'
addpath '../MATLAB functions/time'  

tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
I = imread('Equirectangular.jpg');

%% INITIAL CONDITIONS
date = [2022, 10, 1, 0, 0, 00];
date0 = date2mjd2000(date);
Gr_0 = Get_Greenwich_longitude(date);

mu_earth = astroConstants(13);
R_earth = astroConstants(23);
w_earth = 15.04*pi/(3600*180);

propagation_time = 90*24*3600; % 90 days in seconds

SMA0 = 1.5756e4;    % Km. Given in assignment
ECC0 = 0.5297;      % -. Given in assignment
INC0 = 30.1666;     % degs. Given in assignment
AOP0 = 302.804;    % degs. From Celestrack
RAAN0 = 112.79;    % degs. From Celestrack
TA0 = 162.48414;       % degs. From Celestrack

k = 13;
m = 3;

s0_kep = [SMA0, ECC0, INC0, AOP0, RAAN0, TA0];
T0 = 2*pi*sqrt(SMA0^3/mu_earth);

[r0, v0] = kep2car(s0_kep(1), s0_kep(2), s0_kep(3), s0_kep(5), s0_kep(4), s0_kep(6), mu_earth);
s0_car = [r0', v0'];
%% 1/2.1-INITIAL GROUND TRACK
[t1, s1] = ode113(@(t1, s1) TWOBP(t1, s1, mu_earth), [0, T0], s0_car, options);
[t2, s2] = ode113(@(t2, s2) TWOBP(t2, s2, mu_earth), [0, 24*3600], s0_car, options);
[t3, s3] = ode113(@(t3, s3) TWOBP(t3, s3, mu_earth), [0, 240*3600], s0_car, options);

[~,~, lon1, lat1] =  s2latlong(s1, t1, Gr_0, w_earth);
[~,~, lon2, lat2] =  s2latlong(s2, t2, Gr_0, w_earth);
[~,~, lon3, lat3] =  s2latlong(s3, t3, Gr_0, w_earth);

[lon1, lat1] = ArrangeForPlot(lon1,lat1);
[lon2, lat2] = ArrangeForPlot(lon2,lat2);
[lon3, lat3] = ArrangeForPlot(lon3,lat3);

%% 2.B REPEATING GROUND TRACK
SMA_RGT = GetRGT(w_earth, mu_earth, k, m);
[rRGT, vRGT] = kep2car(SMA_RGT, s0_kep(2), s0_kep(3), s0_kep(5), s0_kep(4), s0_kep(6), mu_earth);
s0_RGT = [rRGT', vRGT'];

[tRGT1, sRGT1] = ode113(@(tRGT, sRGT) TWOBP(tRGT, sRGT, mu_earth), [0, 13*2*pi*sqrt(SMA_RGT^3/mu_earth)], s0_RGT, options);

[~,~, lonRGT1, latRGT1] =  s2latlong(sRGT1, tRGT1, Gr_0, w_earth);

[lonRGT1, latRGT1] = ArrangeForPlot(lonRGT1,latRGT1);

%% 2.C PERTURBATIONS

% For the nominal case only 10 days

[t_perturbed3, s_perturbed3] = ode113(@(t_perturbed, s_perturbed) TwoBodyPerturbed(t_perturbed, s_perturbed, mu_earth, @(t, s) a_J2(t, s) + a_moon(t, s, date0)), [0, 240*3600], s0_car, options);

[~,~, lon_perturbed3, lat_perturbed3] =  s2latlong(s_perturbed3, t_perturbed3, Gr_0, w_earth);

[lon_perturbed3, lat_perturbed3] = ArrangeForPlot(lon_perturbed3,lat_perturbed3);

% For the RGT case only 13 periods

[tRGT_perturbed1, sRGT_perturbed1] = ode113(@(tRGT_perturbed, sRGT_perturbed) TwoBodyPerturbed(tRGT_perturbed, sRGT_perturbed, mu_earth, @(t, s) a_J2(t, s) + a_moon(t, s, date0)), [0, 13*2*pi*sqrt(SMA_RGT^3/mu_earth)], s0_RGT, options);

[~,~, lonRGT_perturbed1, latRGT_perturbed1] =  s2latlong(sRGT_perturbed1, tRGT_perturbed1, Gr_0, w_earth);

[lonRGT_perturbed1, latRGT_perturbed1] = ArrangeForPlot(lonRGT_perturbed1,latRGT_perturbed1);

%% 3 GAUSS EQUATIONS
tic
[t_car_c, s_car_c] = ode113(@(t_perturbed, s_perturbed) TwoBodyPerturbed(t_perturbed, s_perturbed, mu_earth, @(t, s) a_J2(t, s) + a_moon(t,s,date0)), linspace(0, propagation_time, 10001), s0_car, options);
time_car = toc
s0_kep = [s0_kep(1:2),s0_kep(3:6)*pi/180];
tic
[t_kep_c, s_kep_c] = ode113(@(t, kep) Gauss_planetary(t, kep, mu_earth,  @(t, s) a_J2(t, s) + a_moon(t,s,date0)), linspace(0, propagation_time, 10001), s0_kep, options);
time_gauss = toc
car_kep_c = s_car_c*0;
for ii = 1:length(s_kep_c)
    [aa,ee,iii,WW,ww,tt] = ijk2keplerian(s_car_c(ii,1:3)*1000,s_car_c(ii,4:6)*1000);
    if WW > 180
        WW = +WW - 360;
    end
    car_kep_c(ii,:) = [aa/1000,ee,iii*pi/180,ww*pi/180,WW*pi/180,tt*pi/180];
end

%% 6 FILTERING

samples_per_period = length(s_kep_c)*T0/propagation_time;
samples_per_moon_period = length(s_kep_c)*24*28*3600/propagation_time;

filtered_s_kep = Filter(s_kep_c, samples_per_period);                  
filtered_car_kep = Filter(car_kep_c, samples_per_period);
filtered_s_kep_moon = Filter(s_kep_c, samples_per_moon_period);
filtered_car_kep_moon = Filter(car_kep_c, samples_per_moon_period);

%%
figure(31)
hold on
grid on
plot(t_kep_c/(24*3600),s_kep_c(:,1))
plot(t_kep_c/(24*3600),filtered_s_kep(:,1))
plot(t_kep_c/(24*3600),filtered_s_kep_moon(:,1))
title('Short, long and secular evolution of SMA','Interpreter','latex','Fontsize',11)
legend('Short period', 'Long period', 'Secular','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('SMA[KM]','Interpreter','latex','Fontsize',11)

figure(32)
hold on
grid on
plot(t_kep_c/(24*3600),s_kep_c(:,2))
plot(t_kep_c/(24*3600),filtered_s_kep(:,2))
plot(t_kep_c/(24*3600),filtered_s_kep_moon(:,2))
title('Short, long and secular evolution of ECC','Interpreter','latex','Fontsize',11)
legend('Short period', 'Long period', 'Secular','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('Ecc','Interpreter','latex','Fontsize',11)

figure(33)
hold on
grid on
plot(t_kep_c/(24*3600),s_kep_c(:,3)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep(:,3)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep_moon(:,3)*180/pi)
legend('Short period', 'Long period', 'Secular','Interpreter','latex','Fontsize',11)
title('Short, long and secular evolution of INC','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('INC[deg]','Interpreter','latex','Fontsize',11)

figure(34)
hold on
grid on
plot(t_kep_c/(24*3600),s_kep_c(:,4)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep(:,4)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep_moon(:,4)*180/pi)
legend('Short period', 'Long period', 'Secular','Interpreter','latex','Fontsize',11)
title('Short, long and secular evolution of AOP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('AOP[deg]','Interpreter','latex','Fontsize',11)

figure(35)
hold on
grid on
plot(t_kep_c/(24*3600),s_kep_c(:,5)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep(:,5)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep_moon(:,5)*180/pi)
legend('Short period', 'Long period', 'Secular','Interpreter','latex','Fontsize',11)
title('Short, long and secular evolution of RAAN','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('RAAN[deg]','Interpreter','latex','Fontsize',11)

s_kep_c(:,4:6) = mod(s_kep_c(:,4:6),2*pi);       % Arrange the AOP between 0 and 360. In the filtering part it would spoil the average with a jump in the data
filtered_s_kep(:,6) = mod(filtered_s_kep(:,6), 2*pi);
filtered_s_kep_moon(:,6) = mod(filtered_s_kep_moon(:,6), 2*pi);

figure(36)
hold on
grid on
plot(t_kep_c/(24*3600),s_kep_c(:,6)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep(:,6)*180/pi)
plot(t_kep_c/(24*3600),filtered_s_kep_moon(:,6)*180/pi)
legend('Short period', 'Long period', 'Secular','Interpreter','latex','Fontsize',11)
title('Short, long and secular evolution of TA','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('TA[deg]','Interpreter','latex','Fontsize',11)



%% 7 REAL DATA
Horizon_data = importdata("horizons_results_10min.txt");
Horizon_data = Horizon_data.data;
Horizon_kep = zeros(length(Horizon_data),6);
for ii = 1:length(Horizon_data)
    [aa,ee,iii,WW,ww,tt] = ijk2keplerian(Horizon_data(ii,1:3)*1000,Horizon_data(ii,4:6)*1000);
    if WW > 180
        WW = +WW - 360;
    end
    Horizon_kep(ii,:) = [aa/1000,ee,iii*pi/180,ww*pi/180,WW*pi/180,tt*pi/180];
end
Horizon_time = linspace(0, 90*24*3600, 1+90*24*3600/600);

%% PLOTS
%% UNPERTURBED NOMINAL ORBITS
figure(1)
hold on
axis equal
xlabel('longitude [deg]','Interpreter','latex','Fontsize',19)
ylabel('latitude [deg]','Interpreter','latex','Fontsize',19)
xlim([-180, 180])
ylim([-90, 90])
xticks(linspace(-180, 180, 13))
yticks(linspace(-90, 90, 7))
plot(lon1*180/pi, lat1*180/pi, 'r',LineWidth=3)
plot(lon1(1)*180/pi, lat1(1)*180/pi, 'bo',LineWidth=9)
plot(lon1(end)*180/pi, lat1(end)*180/pi, 'bd',LineWidth=9)
h = image(xlim, -ylim, I);
uistack(h, 'bottom')

figure(2)
hold on
axis equal
xlabel('longitude [deg]','Interpreter','latex','Fontsize',19)
ylabel('latitude [deg]','Interpreter','latex','Fontsize',19)
xlim([-180, 180])
ylim([-90, 90])
xticks(linspace(-180, 180, 13))
yticks(linspace(-90, 90, 7))
plot(lon2*180/pi, lat2*180/pi, 'r',LineWidth=3)
plot(lon2(1)*180/pi, lat2(1)*180/pi, 'bo',LineWidth=8)
plot(lon2(end)*180/pi, lat2(end)*180/pi, 'bd',LineWidth=8)
h = image(xlim, -ylim, I);
uistack(h, 'bottom')

figure(3)
hold on
axis equal
xlabel('longitude [deg]','Interpreter','latex','Fontsize',19)
ylabel('latitude [deg]','Interpreter','latex','Fontsize',19)
xlim([-180, 180])
ylim([-90, 90])
xticks(linspace(-180, 180, 13))
yticks(linspace(-90, 90, 7))
plot(lon3*180/pi, lat3*180/pi, 'r',LineWidth=3)
plot(lon3(1)*180/pi, lat3(1)*180/pi, 'bo',LineWidth=8)
plot(lon3(end)*180/pi, lat3(end)*180/pi, 'bd',LineWidth=8)
h = image(xlim, -ylim, I);
uistack(h, 'bottom')

%% PERTURBED NOMINAL 

figure(4)
hold on
axis equal
xlabel('longitude [deg]','Interpreter','latex','Fontsize',19)
ylabel('latitude [deg]','Interpreter','latex','Fontsize',19)
xlim([-180, 180])
ylim([-90, 90])
xticks(linspace(-180, 180, 13))
yticks(linspace(-90, 90, 7))
plot(lon_perturbed3*180/pi, lat_perturbed3*180/pi, 'r',LineWidth=3)
plot(lon_perturbed3(1)*180/pi, lat_perturbed3(1)*180/pi, 'bo',LineWidth=8)
plot(lon_perturbed3(end)*180/pi, lat_perturbed3(end)*180/pi, 'bd',LineWidth=8)
h = image(xlim, -ylim, I);
uistack(h, 'bottom')

%% UNPERTURBED RGT

figure(5)
hold on
axis equal
xlabel('longitude [deg]','Interpreter','latex','Fontsize',19)
ylabel('latitude [deg]','Interpreter','latex','Fontsize',19)
xlim([-180, 180])
ylim([-90, 90])
xticks(linspace(-180, 180, 13))
yticks(linspace(-90, 90, 7))
plot(lonRGT1*180/pi, latRGT1*180/pi, 'r',LineWidth=3)
plot(lonRGT1(1)*180/pi, latRGT1(1)*180/pi, 'bo',LineWidth=8)
plot(lonRGT1(end)*180/pi, latRGT1(end)*180/pi, 'bd',LineWidth=8)
h = image(xlim, -ylim, I);
uistack(h, 'bottom')

%% PERTURBED RGT

figure(6)
hold on
axis equal
xlabel('longitude [deg]','Interpreter','latex','Fontsize',19)
ylabel('latitude [deg]','Interpreter','latex','Fontsize',19)
xlim([-180, 180])
ylim([-90, 90])
xticks(linspace(-180, 180, 13))
yticks(linspace(-90, 90, 7))
plot(lonRGT_perturbed1*180/pi, latRGT_perturbed1*180/pi, 'r',LineWidth=3)
plot(lonRGT_perturbed1(1)*180/pi, latRGT_perturbed1(1)*180/pi, 'bo',LineWidth=8)
plot(lonRGT_perturbed1(end)*180/pi, latRGT_perturbed1(end)*180/pi, 'bd',LineWidth=8)
h = image(xlim, -ylim, I);
uistack(h, 'bottom')

%% Other plots

figure(61)
hold on
grid on
title('SMA evolution simulated vs real data','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,1), 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,1),'r-')
plot(Horizon_time/(24*3600), Horizon_kep(:,1), 'm')
legend('Gauss equations', 'Perturbed 2BP', 'Real data','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('SMA[KM]','Interpreter','latex','Fontsize',11)

figure(62)
hold on
grid on
title('ECC evolution simulated vs real data','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,2), 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,2), 'r-')
plot(Horizon_time/(24*3600), Horizon_kep(:,2), 'm')
legend('Gauss equations', 'Perturbed 2BP', 'Real data','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('Ecc','Interpreter','latex','Fontsize',11)

figure(63)
hold on
grid on
title('INC evolution simulated vs real data','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,3)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,3)*180/pi, 'r-')
plot(Horizon_time/(24*3600), Horizon_kep(:,3)*180/pi, 'm')
legend('Gauss equations', 'Perturbed 2BP', 'Real data')
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('INC[deg]','Interpreter','latex','Fontsize',11)

figure(64)
hold on
grid on
title('AOP evolution simulated vs real data','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,4)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,4)*180/pi, 'r-')
plot(Horizon_time/(24*3600), Horizon_kep(:,4)*180/pi, 'm')
legend('Gauss equations', 'Perturbed 2BP', 'Real data','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('AOP[deg]','Interpreter','latex','Fontsize',11)

figure(65)
hold on
grid on
title('RAAN evolution simulated vs real data','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,5)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,5)*180/pi, 'r-')
plot(Horizon_time/(24*3600), Horizon_kep(:,5)*180/pi, 'm')
legend('Gauss equations', 'Perturbed 2BP', 'Real data','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('RAAN[deg]','Interpreter','latex','Fontsize',11)

figure(66)
hold on
grid on
title('TA evolution simulated vs real data','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,6)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,6)*180/pi, 'r-')
plot(Horizon_time/(24*3600), Horizon_kep(:,6)*180/pi, 'm')
legend('Gauss equations', 'Perturbed 2BP', 'Real data','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('TA[deg]','Interpreter','latex','Fontsize',11)


figure(71)
hold on
grid on
title('SMA according to Gauss equations and TBP propagation','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,1), 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,1),'r-')
legend('Gauss equations', 'Perturbed 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('SMA[KM]','Interpreter','latex','Fontsize',11)

figure(72)
hold on
grid on
title('ECC according to Gauss equations and TBP propagation','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,2), 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,2), 'r-')
legend('Gauss equations', 'Perturbed 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('Ecc','Interpreter','latex','Fontsize',11)

figure(73)
hold on
grid on
title('INC according to Gauss equations and TBP propagation','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,3)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,3)*180/pi, 'r-')
legend('Gauss equations', 'Perturbed 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('INC[deg]','Interpreter','latex','Fontsize',11)

figure(74)
hold on
grid on
title('AOP according to Gauss equations and TBP propagation','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,4)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,4)*180/pi, 'r-')
legend('Gauss equations', 'Perturbed 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('AOP[deg]','Interpreter','latex','Fontsize',11)

figure(75)
hold on
grid on
title('RAAN according to Gauss equations and TBP propagation','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,5)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,5)*180/pi, 'r-')
legend('Gauss equations', 'Perturbed 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('RAAN[deg]','Interpreter','latex','Fontsize',11)

figure(76)
hold on
grid on
title('TA according to Gauss equations and TBP propagation','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), s_kep_c(:,6)*180/pi, 'b-')
plot(t_car_c/(24*3600), car_kep_c(:,6)*180/pi, 'r-')
legend('Gauss equations', 'Perturbed 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('TA[deg]','Interpreter','latex','Fontsize',11)

figure(81)
hold on
grid on
title('Error in SMA between Gauss and 2BP','Interpreter','latex','Fontsize',11)
plot(t_kep_c/(24*3600), (s_kep_c(:,1)-car_kep_c(:,1))/SMA0, 'b-')
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('$(SMA_g - SMA_t)/SMA_0$','Interpreter','latex','Fontsize',11)

figure(82)
hold on
grid on
plot(t_kep_c/(24*3600), s_kep_c(:,2)-car_kep_c(:,2), 'b-')
title('Error in ECC between Gauss and 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('$ECC_g - ECC_t$','Interpreter','latex','Fontsize',11)

figure(83)
hold on
grid on
plot(t_kep_c/(24*3600), (s_kep_c(:,3)*180/pi-car_kep_c(:,3)*180/pi)/360, 'b-')
title('Error in INC between Gauss and 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('$(INC_g - INC_t)/360$','Interpreter','latex','Fontsize',11)

figure(84)
hold on
grid on
plot(t_kep_c/(24*3600), (s_kep_c(:,4)*180/pi-car_kep_c(:,4)*180/pi)/360, 'b-')
title('Error in AOP between Gauss and 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('$(AOP_g - AOP_t)/360$','Interpreter','latex','Fontsize',11)

figure(85)
hold on
grid on
plot(t_kep_c/(24*3600), (s_kep_c(:,5)*180/pi-car_kep_c(:,5)*180/pi)/360, 'b-')
title('Error in RAAN between Gauss and 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('$(RAAN_g - RAAN_t)/360$','Interpreter','latex','Fontsize',11)

figure(86)
hold on
grid on
plot(t_kep_c/(24*3600), (s_kep_c(:,6)*180/pi-car_kep_c(:,6)*180/pi)/360, 'b-')
title('Error in TA between Gauss and 2BP','Interpreter','latex','Fontsize',11)
xlabel('Elapsed days','Interpreter','latex','Fontsize',11)
ylabel('$(TA_g - TA_t)/360$','Interpreter','latex','Fontsize',11)

tcolor = floor(t_car_c/T0);
tcolor(end+1) = NaN;
s_car_c(end+1,:) = NaN;

colormap('jet')

figure (100)
hold on
grid on
axis equal
xlabel('x (Km)','Interpreter','latex','Fontsize',11)
ylabel('y (Km)','Interpreter','latex','Fontsize',11)
zlabel('z (Km)','Interpreter','latex','Fontsize',11)
title('Orbit representation','Interpreter','latex','Fontsize',11)
DrawEarth
patch(s_car_c(:,1), s_car_c(:,2), s_car_c(:,3), tcolor, 'EdgeColor', 'interp')
hcb = colorbar;
title(hcb, 'Orbit number','Interpreter','latex','Fontsize',11)

s_car_c(end,:) = [];
