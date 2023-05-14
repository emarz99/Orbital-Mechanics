%% Orbital Mechanics - Assignment 1 - Interplanetary explorer mission
clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% Introduction

% Assignement 1 Data
dep_Planet = 1; % Mercury
fb_Planet = 2;  % Venus
arr_Planet = 3; % Earth

dep_mu = astroConstants(10+dep_Planet);
fb_mu = astroConstants(10+fb_Planet);
arr_mu = astroConstants(10+arr_Planet);
mu_Sun = astroConstants(4);

earliest_dep_date = [2027 5 1 0 0 0];
latest_arr_date = [2067 5 1 0 0 0];

earliestDepMJD2000 = date2mjd2000(earliest_dep_date);
latestArrMJD2000 = date2mjd2000(latest_arr_date);


[kep_Me,~] = uplanet(date2mjd2000(earliest_dep_date), dep_Planet);
[kep_Ve,~] = uplanet(date2mjd2000(earliest_dep_date), fb_Planet);
[kep_Ea,~] = uplanet(date2mjd2000(earliest_dep_date), arr_Planet);

T_Me = 2*pi*sqrt(kep_Me(1)^3/mu_Sun);
T_Ve = 2*pi*sqrt(kep_Ve(1)^3/mu_Sun);
T_Ea = 2*pi*sqrt(kep_Ea(1)^3/mu_Sun);
Tsy_MeVe = ((1/T_Me - 1/T_Ve)^-1)/(3600*24);
Tsy_VeEa = ((1/T_Ve - 1/T_Ea)^-1)/(3600*24);

%   Mercury and Venus will approximately return to theirs respective positions every 5.5 years / 14 synodic periods
%   Venus and Earth will return to theirs respective positions every 8 years / 5 synodic periods

%% Lambert tranfers

N  = 1462;       % every 10 days (14611 for each day)
MJD = linspace(earliestDepMJD2000,latestArrMJD2000,N);

[min_DV1, ToF_minDV1, depdate_minDV1, arrdate_minDV1, ~, all_DV1, ~, ~, v1_fb] = lambertSolver(earliest_dep_date,latest_arr_date, earliest_dep_date, latest_arr_date, dep_Planet, fb_Planet, mu_Sun, N, 1);
[min_DV2, ToF_minDV2, depdate_minDV2, arrdate_minDV2, ~, ~, all_DV2, v2_fb, ~] = lambertSolver(earliest_dep_date,latest_arr_date, earliest_dep_date, latest_arr_date, fb_Planet, arr_Planet, mu_Sun, N, 2);

%% PorkChop Plot
toll_1 = 0.65;
orbT1=7;          %[years] seems to be the repetitionPeriod of firstTransf by graphic approach
toll_2 = 0.5;
orbT2=10;         %[years] seems to be the repetitionPeriod of the secondTransf by graphic approach

PorkChop(earliest_dep_date,latest_arr_date, all_DV1, N, toll_1, orbT1);
PorkChop(earliest_dep_date,latest_arr_date, all_DV2, N, toll_2, orbT2);

%% min(DV1+DV2) Solution

for i = 1:N
    for j = i+1:N
        for k = j+1:N
            if ~isnan(all_DV1(i,j)) && ~isnan(all_DV2(j,k)) && (abs(all_DV1(i,j) - min_DV1)+abs(all_DV2(j,k) - min_DV2)) < 2

                dep = mjd20002date(MJD(i));
                fb  = mjd20002date(MJD(j));
                arr = mjd20002date(MJD(k));
                [DV_fb,~,~,~] = flyby(v1_fb(i,:,j) , v2_fb(j,:,k), fb_Planet, MJD(j));

                if isnan(DV_fb)
                    warning('min(DV1+DV2) soultion impossible due to gravity assist maneuver too close to planet');
                else
                    disp('Solution found');
                end

                return
            end
        end
    end
end

%% min(DV1+DV2) finite solution

for toll = 1:10
    for i = 1:N
        for j = i+1:N
            for k = j+1:N
                if ~isnan(all_DV1(i,j)) && ~isnan(all_DV2(j,k)) && (abs(all_DV1(i,j) - min_DV1)+abs(all_DV2(j,k) - min_DV2)) < toll && (abs(all_DV1(i,j) - min_DV1)+abs(all_DV2(j,k) - min_DV2)) >= toll-1

                    dep = mjd20002date(MJD(i));
                    fb  = mjd20002date(MJD(j));
                    arr = mjd20002date(MJD(k));
                    [DV_fb,~,~,~] = flyby( v1_fb(i,:,j), v2_fb(j,:,k), fb_Planet, MJD(j));

                    if ~isnan(DV_fb)

                        disp('Acceptable min(DV1+DV2) soultion found')
                        Dv_tot = all_DV1(i,j) + all_DV2(j,k) + DV_fb;  
                        DV1 = all_DV1(i,j);                            
                        DV2 = all_DV2(j,k);                            
                        return

                    end

                end
            end
        end
    end
end

%% Stem plots

DV_tot = NaN(N,N,N);
DV_fb = NaN(N,N,N);
DV1_DV2 = NaN(N,N,N);

for i = 1:N
    for j = i+1:N
        for k = j+1:N
            if ~isnan(all_DV1(i,j)) && ~isnan(all_DV2(j,k))
            
                [DV_fb(i,j,k), ~, ~, ~] = flyby( v1_fb(i,:,j), v2_fb(j,:,k), fb_Planet, MJD(j));
                DV1_DV2(i,j,k) = all_DV1(i,j) + all_DV2(j,k);
                DV_tot(i,j,k) = DV1_DV2(i,j,k) + DV_fb(i,j,k);
                
            end
        end
    end
end

figure()
hold on
for j = 1:N
    for k = 1:N
    stem(DV_tot(:,j,k), DV1_DV2(:,j,k), 'b')
    end
end
xlabel('Total cost of the mission [Km/s]')
ylabel('Departure plus arrival cost [Km/s]')
title('Costs comparsion')

figure()
hold on
for j = 1:N
    for k = 1:N
    stem(DV_tot(:,j,k), DV_fb(:,j,k), 'r')
    end
end
xlabel('Total cost of the mission [Km/s]')
ylabel('Fly-by cost [Km/s]')
title('Costs comparsion')

figure()
hold on
xlabel('Total cost of the mission [Km/s]')
ylabel('ToFs [days]')
title('ToF / DVtot comparsion')
for i = 1:N
        for k = i+1:N
            p1=stem(DV_tot(i,:,k), MJD(:)-MJD(i), 'b');
            p2=stem(DV_tot(i,:,k), -(MJD(k)-MJD(:)), 'r');
        end
end
xlim([0 30]);
legend([p1 p2],'Mercury-Venus ToF','Venus-Earth ToF')


%% Grid search
DV_tot = NaN(N,N);
fb_MJD = NaN(N,N);

for i = 1:N
    for j = i+1:N
        for k = j+1:N

            if ~isnan(all_DV1(i,j)) && ~isnan(all_DV2(j,k))
                [DV_fb, ~, ~, ~] = flyby( v1_fb(i,:,j), v2_fb(j,:,k), fb_Planet, MJD(j));
                DVtot = all_DV1(i,j) + all_DV2(j,k) + DV_fb;

                if ~isnan(DVtot) && isnan(DV_tot(i,k))
                    fb_MJD(i,k) = MJD(j);
                    DV_tot(i,k) = DVtot;
                end

                if ~isnan(DVtot) && ~isnan(DV_tot(i,k)) && DVtot < DV_tot(i,k)

                    fb_MJD(i,k) = MJD(j);
                    DV_tot(i,k) = DVtot;
                end
            end

        end
    end
end

min_DVtot_guess = min(DV_tot,[],'all')
[idx1, idx2] = find(DV_tot==min_DVtot_guess);


toll = 0.75;
PorkChop(earliest_dep_date,latest_arr_date, DV_tot, N, toll);


earliest_dep_date = mjd20002date(MJD(idx1-1));
latest_dep_date   = mjd20002date(MJD(idx1+1));

earliest_fb_date = mjd20002date(fb_MJD(idx1,idx2)-(MJD(2)-MJD(1)));
latest_fb_date   = mjd20002date(fb_MJD(idx1,idx2)+(MJD(2)-MJD(1)));

earliest_arr_date = mjd20002date(MJD(idx2-1));
latest_arr_date   = mjd20002date(MJD(idx2+1));


%% Refined search
N = 20*24;    % every hour (for grid search every 10) 

MJD1 = linspace(date2mjd2000(earliest_dep_date), date2mjd2000(latest_dep_date), N);
MJD2 = linspace(date2mjd2000(earliest_fb_date), date2mjd2000(latest_fb_date), N);
MJD3 = linspace(date2mjd2000(earliest_arr_date), date2mjd2000(latest_arr_date), N);

[min_DV1, ToF_minDV1, depdate_minDV1, arrdate_minDV1, ~, all_DV1, ~, v1_t1, v1_fb] = lambertSolver(earliest_dep_date, latest_dep_date, earliest_fb_date, latest_fb_date, dep_Planet, fb_Planet, mu_Sun, N, 1);
[min_DV2, ToF_minDV2, depdate_minDV2, arrdate_minDV2, ~, ~, all_DV2, v2_fb, ~] = lambertSolver(earliest_fb_date, latest_fb_date, earliest_arr_date, latest_arr_date, fb_Planet, arr_Planet, mu_Sun, N, 2);

DV_tot = NaN(N,N,N);
DV_fb = NaN(N,N,N);
DV1_DV2 = NaN(N,N,N);

for i = 1:N
    for j = 1:N
        for k = 1:N
            if ~isnan(all_DV1(i,j)) && ~isnan(all_DV2(j,k))

                [DV_fb(i,j,k), ~, ~, ~] = flyby( v1_fb(i,:,j), v2_fb(j,:,k), fb_Planet, MJD2(j));
                DV1_DV2(i,j,k) = all_DV1(i,j) + all_DV2(j,k);
                DV_tot(i,j,k) = DV1_DV2(i,j,k) + DV_fb(i,j,k);

            end
        end
    end
end

min_DVtot = min(DV_tot,[],'all');
[idx1, idx2, idx3] = ind2sub(size(DV_tot), find(DV_tot==min_DVtot));

DV1 = all_DV1(idx1,idx2);
DV2 = all_DV2(idx2,idx3);
DVfb = DV_fb(idx1,idx2,idx3);

v1_arr = v1_t1(idx1,:,idx2);

depMJD_minDV = MJD1(idx1);
fbMJD_minDV  = MJD2(idx2);
arrMJD_minDV = MJD3(idx3);

depdate_minDV = mjd20002date(depMJD_minDV);
fbdate_minDV  = mjd20002date(fbMJD_minDV);
arrdate_minDV = mjd20002date(arrMJD_minDV);


%% Mission Plots and Cost Definition

% depdate_minDV = [2050 2 18 16 0 0];
% fbdate_minDV  = [2050 6 6 2 0 0];
% arrdate_minDV = [2050 11 25 12 0 0];



[v_mer, vdep_t, r_mer, varr_t, v_earth, r_earth,  v_minus, v_plus, r_venus, v_inf_i, v_inf_f, Dv_ga_vec, Time_in_SOI, hp, delta_deg]=mission_plot(dep_Planet, depdate_minDV, fb_Planet, fbdate_minDV, arr_Planet, arrdate_minDV, 1);

DV1_vec=vdep_t-v_mer
DV2_vec=v_earth-varr_t
DV_fb_vec=v_plus-v_minus
DV_ga_vec=Dv_ga_vec

DV1= norm(DV1_vec)
DV2=norm(DV2_vec)
DV_effective=norm(DV_fb_vec)           %==DV_effective
DV_cost_peric=norm(DV_ga_vec)          %==DV_cost_peric

DV_tot=DV1+DV2+DV_cost_peric;

%% Double flyby

N = round(4*T_Me/(24*3600));                % 4 Me revolutions

MJD = (depMJD_minDV-N+1:1:depMJD_minDV);

r_dep = zeros(N,3);
v_dep = zeros(N,3);

for i = 1:N

    [kep_Planet,~] = uplanet(MJD(i),dep_Planet);
    [r_dep(i,:), v_dep(i,:)] = kep2car(kep_Planet(1),kep_Planet(2),kep_Planet(3),kep_Planet(4),kep_Planet(5),kep_Planet(6),mu_Sun);

end

r_arr = r_dep(end,:);

for j = 0:4
    for k = 0:1

        orbitType = 0;
        Nrev = j;
        Ncase = k;
        optionsLMR = 1;

        partial_DVtot = NaN(N,1);
        DV_fb = NaN(N,1);
        DV_1 = NaN(N,1);
        v1_t = NaN(N,3);
        v2_t = NaN(N,3);

        for i = 1:N-1

            ToF = (MJD(end)-MJD(i))*24*3600;
            [~,~,~,~,v1_t(i,:),v2_t(i,:),~,~] = lambertMR(r_dep(i,:),r_arr,ToF,mu_Sun,orbitType,Nrev,Ncase,optionsLMR);
            DV_1(i) = norm(v1_t(i,:) - v_dep(i,:));
            [DV_fb(i), ~, ~, ~] = flyby( v2_t(i,:), v1_arr, dep_Planet, MJD(end));

            if DV_1(i) > 0.5
                partial_DVtot(i) = DV_1(i) + DV_fb(i);
            end

        end

        if ~isnan(min(partial_DVtot))

            idx = find(partial_DVtot==min(partial_DVtot));
            newDV1 = DV_1(idx);
            DV_fb1 = DV_fb(idx,:);
            depMJD_new = MJD(idx);

            new_depdate = mjd20002date(depMJD_new);
            min_DVtot = newDV1 + DV_fb1 + DVfb + DV2;

            mission_plot(dep_Planet, depdate_minDV, fb_Planet, fbdate_minDV, arr_Planet, arrdate_minDV, 1, new_depdate,[Nrev Ncase])
            return

        end


    end
end