function [deltaV_total,deltaV1,V_initial,deltaV2,V_final] = delta_V_Interplanetary(t1,t2,Celestial_1,Celestial_2)
% This function calculates the deltaV of an interplanetary heliocentric
% mission using the Lambert Solver between Planet 1 and Planet 2 given by 
% their indeces, as seen in the uplanet and ephNEO function, for a specific
% departure and arrival date.
%
% PROTOTYPE
%     [deltaV,deltaV1,VI,deltaV2,VF] = delta_V_Interplanetary(t1,t2,ID_1,ID_2)
% 
% INPUT:
%     t1[1]             Departure Time given in mjd2000 format, as seen in
%                       the date2mjd2000 function       [T]
%     t2[1]             Arrival Time given in mjd2000 format, as seen in 
%                       the date2mjd2000 function       [T]
%     Celestial_1[1]    ID of Departure Planet as described in the uplanet
%                       function                        [-]
%     Celestial_2[1]    ID of Arrival Planet as described in the uplanet 
%                       function                        [-]
%
% OUTPUT:
%     deltaV_total[1]   deltaV output of the mission given as the norm of
%                       the deltaV1 and deltaV2                 [L/T]
%     deltaV1[3]        deltaV vector of the first manoeuvre    [L/T]
%     V_initial[3]      Velocity at t1                          [L/T]
%     deltaV2[3]        deltaV vector of the second manoeuvre   [L/T]
%     V_final[3]        Velocity at t2                          [L/T]
%
% CONTRIBUTORS:
%     Carlos Albiñana
%
% VERSIONS:
%     2022-12-9:    First Version
%     2022-12-10:   Added deltaV1 and deltaV2 vectors as outputs
%     2022-12-12:   Debugging for Missmatch in Celestial IDs   
%     2022-12-14:   Added V_initial and V_final outputs
%
% CALLED FUNCTIONS:
%     WEBEEP:
%     uplanet
%     ephNEO
%     lambertMR
%     
%     OWN FUNCTIONS:
%     kep2car
% IDs:
% uplanet function
% Integer number identifying the celestial body (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
% ID>11 is used for ephNEO function


%% Position and Velocity of Planets at t1 and t2

% uplanet function
% Keplerian Planets
%  OUTPUT:
%	kep[6]    	Mean Keplerian elements of date
%                 kep = [a e i Om om theta] [km, rad]
%	ksun[1]     Gravity constant of the Sun [km^3/s^2]

%% Gravitational constant of the Sun [km^3/s^2]
G=6.67259e-20;
msun=1.988919445342813e+030;
ksun=msun*G;

%% Position and Velocity of Planets at t1

% Keplerian elements    - Planet 1 - Arrival Time
% Position and Velocity - Planet 1 - Arrival Time

if Celestial_1<11
    [Celestial_1_kep_t1,~]=uplanet(t1,Celestial_1);
else
    [Celestial_1_kep_t1,~,~,~]=ephNEO(t1,Celestial_1);
end

[ys_Celestial1_t1,v_Celestial1_t1]=kep2car(Celestial_1_kep_t1(1),...
    Celestial_1_kep_t1(2),...
    rad2deg(Celestial_1_kep_t1(3)),...
    rad2deg(Celestial_1_kep_t1(4)),...
    rad2deg(Celestial_1_kep_t1(5)),...
    rad2deg(Celestial_1_kep_t1(6)),...
    ksun);

%% Position and Velocity of Planets at t2

% Keplerian elements    - Planet 2 - Arrival Time
% Position and Velocity - Planet 2 - Arrival Time 

if Celestial_2<11
    [Celestial_2_kep_t2,~]=uplanet(t2,Celestial_2);
else
    [Celestial_2_kep_t2,~,~,~]=ephNEO(t2,Celestial_2);
end

[ys_Ceslestial2_t2,v_Celestial2_t2]=kep2car(Celestial_2_kep_t2(1),...
    Celestial_2_kep_t2(2),...
    rad2deg(Celestial_2_kep_t2(3)),...
    rad2deg(Celestial_2_kep_t2(4)),...
    rad2deg(Celestial_2_kep_t2(5)),...
    rad2deg(Celestial_2_kep_t2(6)),...
    ksun);


%% Lambert solution 

ToF=(t2-t1)*24*3600; % CONVERT INTO SECONDS. DAYS TO SECONDS

% LAMBERT SOLVER
[~,~,~,~,V_initial,V_final,~,~] = lambertMR( ys_Celestial1_t1, ys_Ceslestial2_t2, ToF, ksun, 0, 0, 0 );


deltaV1=V_initial-v_Celestial1_t1;
deltaV2=V_final-v_Celestial2_t2;

deltaV_total=vecnorm(deltaV1)+vecnorm(deltaV2);


end