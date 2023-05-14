function [minDV, ToF_minDV, depdate_minDV, arrdate_minDV, DVtot, DV_1, DV_2, v1_t, v2_t] = lambertSolver(earliest_dep_date,latest_dep_date, earliest_arr_date, latest_arr_date, dep_Planet, arr_Planet, mu_Sun, N, option)


% mission_plot.m - collect all the feasibles mission costs
%
% DESCRIPTION:taking as input arrival and departure windows built matrixes
%            containing numerical values of all feasible transfers 

%INPUT:
%earliest_dep_date,latest_dep_date[6]     dates composing departure window
%earliest_arr_date, latest_arr_date[6]    dates composing arrival window
% dep_Planet, arr_Planet [1]   Integer number identifying the body 
%                                   1:   Mercury
%                                   2:   Venus
%                                   3:   Earth
%                                   4:   Mars
%                                   5:   Jupiter
%                                   6:   Saturn
%                                   7:   Uranus
%                                   8:   Neptune
%                                   9:   Pluto
%                                   10:  Sun
% mu_Sun,                   astronomic constant of Sun [km^3/s^2]
% N                                 number of points considered 
% option
% Contour Plot options:     = 0   :   min DV= min(DVtot)   
%                           = 1   :   min DV= min(DV1)
%                           = 2   :   min DV= min(DV2)

%OUTPUT:
% minDV          [1]         min value of DV ---> funcion of the selected option!
% ToF_minDV      [1]         [days] time of flight corresponding to min DV
% depdate_minDV  [1x6]       Gregorian date of departure associated to minDV
% arrdate_minDV  [1x6]       Gegorian date of departure associated to minDV
% DVtot          [NxN]       Matrix containing all the feasible transfer cost
% DV_1           [NxN]       Matrix containing all the 1st maneuver costs
% DV_2           [NxN]       Matrix containing all the 1st maneuver costs
% v1_t           [Nx3xN]     3D matrix containing vector of the s/c
%                            velocity on the first point of the transfer arc
% v2_t           [Nx3xN]    3D matrix containing vector of the s/c
%                           velocity on the last point of the transfer arc

% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

% Departure and arrival time windows
initialDepMJD2000 = date2mjd2000(earliest_dep_date);
initialArrMJD2000 = date2mjd2000(earliest_arr_date);

finalDepMJD2000 = date2mjd2000(latest_dep_date);
finalArrMJD2000 = date2mjd2000(latest_arr_date);

depMJD = linspace(initialDepMJD2000,finalDepMJD2000,N);
arrMJD = linspace(initialArrMJD2000,finalArrMJD2000,N);

% Planets' orbits
r_dep_Planet = zeros(N,3);
v_dep_Planet = zeros(N,3);
r_arr_Planet = zeros(N,3);
v_arr_Planet = zeros(N,3);

% Lambert Option 
orbitType = 0;       
Nrev = 0;
Ncase = 0;
optionsLMR = 1;

% Lambert output
DVtot = NaN(N,N);
DV_1 = NaN(N,N);
DV_2 = NaN(N,N);
v1_t = NaN(N,3,N);
v2_t = NaN(N,3,N);

for k = 1:N
    [kep_dep_Planet,~] = uplanet(depMJD(k),dep_Planet);
    [r_dep_Planet(k,:),v_dep_Planet(k,:)] = kep2car(kep_dep_Planet(1),kep_dep_Planet(2),kep_dep_Planet(3),kep_dep_Planet(4),kep_dep_Planet(5),kep_dep_Planet(6),mu_Sun);
    
    [kep_arr_Planet,~] = uplanet(arrMJD(k),arr_Planet);
    [r_arr_Planet(k,:),v_arr_Planet(k,:)] = kep2car(kep_arr_Planet(1),kep_arr_Planet(2),kep_arr_Planet(3),kep_arr_Planet(4),kep_arr_Planet(5),kep_arr_Planet(6),mu_Sun);
end


for j = 1:N
    for k = 1:N
        if arrMJD(k)-depMJD(j) <= 365 && arrMJD(k)-depMJD(j) > 0
            
            ToF = (arrMJD(k)-depMJD(j))*24*3600;
            [~,~,~,~,v1_t(j,:,k),v2_t(j,:,k),~,~] = lambertMR(r_dep_Planet(j,:),r_arr_Planet(k,:),ToF,mu_Sun,orbitType,Nrev,Ncase,optionsLMR);
            DV_1(j,k) = norm(v1_t(j,:,k) - v_dep_Planet(j,:));
            DV_2(j,k) = norm(v_arr_Planet(k,:) - v2_t(j,:,k));
            DVtot(j,k) = DV_1(j,k) + DV_2(j,k);
             
        end        
    end
end

switch option
    case 0        
        DV_plot = DVtot;
    case 1
        DV_plot = DV_1; 
    case 2
        DV_plot = DV_2;         
end

minDV = min(DV_plot,[],'all');
[idx1, idx2] = find(DV_plot == minDV);

depMJD_minDV = depMJD(idx1);
arrMJD_minDV = arrMJD(idx2);

ToF_minDV = arrMJD_minDV-depMJD_minDV;

depdate_minDV = mjd20002date(depMJD_minDV);
arrdate_minDV = mjd20002date(arrMJD_minDV);

return