function[v_dep_planet, vdep_t, r_dep_Planet, varr_t, v_arr_planet, r_arr_Planet,  V_minus, V_plus, r_fb_Planet, v_inf_i, v_inf_f, DV_ga_vec, Time_in_SOI, hp, delta_deg]=mission_plot(dep_Planet, depDate, fb_Planet, fbDate, arr_Planet, arrDate, option, ga2Date, ga2option)


% mission_plot.m - collect all the mission plots
%
% DESCRIPTION:taking as input the choosen dates of the mission plots the
% transfer arcs, the fly-bys and the definition of the maneuver points

%INPUT:
% dep_Planet , fb_Planet,  arr_Planet   [1]   Integer number identifying the bodies 
%                                                1:   Mercury
%                                                2:   Venus
%                                                3:   Earth
%                                                4:   Mars
%                                                5:   Jupiter
%                                                6:   Saturn
%                                                7:   Uranus
%                                                8:   Neptune
%                                                9:   Pluto
%                                                10:  Sun   
% depDate, fbDate, arrDate           [1x6]   Date in the Gregorian calendar, as a 6-element vector
%                                              [year, month, day, hour, minute, second]  
% options:    = 1 :  Static plot
%             = 2 :  Dynamic plot
%optional INPUT:
%    ga2Date                         [1x6] Date in the Gregorian calendar, as a 6-element vector of 2Fly-by date
%                                              [year, month, day, hour, minute, second]
                             
%OUTPUT:
%v_dep_planet      [3]    vel before 1st maneuver
% vdep_t           [3]    vel after 1st maneuver
% r_dep_Planet     [3]    position of 1st meneuver
% varr_t           [3]    velocity before last maneuver  
% v_arr_planet     [3]    velocity after last maneuver
% r_arr_Planet     [3]    position of last maneuver
% V_minus          [3]    velocity before fly-by maneuver
% V_plus           [3]    velocity after fly-by maneuver
% r_fb_Planet      [3]    position of fly-by maneuver
% v_inf_i          [3]    velocity before fly-by maneuver in venus ref frame
% v_inf_f          [3]    velocity after fly-by maneuver in venus ref frame
% DV_ga_vec        [3]    cost of fly-by (tangentDV at hyperbola pericenter)
% Time_in_SOI      [1]    [s] time spent in Venus SOI
% hp               [1]    pericenter altitude of hyperbola wrt venus
% delta_deg        [1]    angle between v_inf_i and v_inf_f


% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

mu_Sun = astroConstants(4);
AU = astroConstants(2);

% Lambert options
orbitType = 0;       
Nrev = 0;
Ncase = 0;
options = 1;
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-13);


depTime = date2mjd2000(depDate);
fbTime  = date2mjd2000(fbDate);
arrTime = date2mjd2000(arrDate);

kep_dep = uplanet( depTime, dep_Planet);
kep_fb  = uplanet( fbTime, fb_Planet);
kep_arr = uplanet( arrTime, arr_Planet);

[r_dep_Planet, v_dep_Planet] = kep2car(kep_dep(1),kep_dep(2),kep_dep(3),kep_dep(4),kep_dep(5),kep_dep(6),mu_Sun);
[r_fb_Planet, ~] = kep2car(kep_fb(1),kep_fb(2),kep_fb(3),kep_fb(4),kep_fb(5),kep_fb(6),mu_Sun);
[r_arr_Planet, v_arr_Planet] = kep2car(kep_arr(1),kep_arr(2),kep_arr(3),kep_arr(4),kep_arr(5),kep_arr(6),mu_Sun);

v_dep_planet=v_dep_Planet';
v_arr_planet=v_arr_Planet';

kep_dep_dep = uplanet( depTime, dep_Planet);
kep_dep_fb  = uplanet( fbTime, dep_Planet);
kep_dep_arr = uplanet( arrTime, dep_Planet);
[r_dep_dep, ~] = kep2car(kep_dep_dep(1),kep_dep_dep(2),kep_dep_dep(3),kep_dep_dep(4),kep_dep_dep(5),kep_dep_dep(6),mu_Sun);
[r_dep_fb, ~] = kep2car(kep_dep_fb(1),kep_dep_fb(2),kep_dep_fb(3),kep_dep_fb(4),kep_dep_fb(5),kep_dep_fb(6),mu_Sun);
[r_dep_arr, ~] = kep2car(kep_dep_arr(1),kep_dep_arr(2),kep_dep_arr(3),kep_dep_arr(4),kep_dep_arr(5),kep_dep_arr(6),mu_Sun);


kep_fb_dep = uplanet( depTime, fb_Planet);
kep_fb_fb  = uplanet( fbTime, fb_Planet);
kep_fb_arr = uplanet( arrTime, fb_Planet);
[r_fb_dep, v_fb_dep] = kep2car(kep_fb_dep(1),kep_fb_dep(2),kep_fb_dep(3),kep_fb_dep(4),kep_fb_dep(5),kep_fb_dep(6),mu_Sun);
[r_fb_fb, ~] = kep2car(kep_fb_fb(1),kep_fb_fb(2),kep_fb_fb(3),kep_fb_fb(4),kep_fb_fb(5),kep_fb_fb(6),mu_Sun);
[r_fb_arr, ~] = kep2car(kep_fb_arr(1),kep_fb_arr(2),kep_fb_arr(3),kep_fb_arr(4),kep_fb_arr(5),kep_fb_arr(6),mu_Sun);


kep_arr_dep = uplanet( depTime, arr_Planet);
kep_arr_fb  = uplanet( fbTime, arr_Planet);
kep_arr_arr = uplanet( arrTime, arr_Planet);
[r_arr_dep, v_arr_dep] = kep2car(kep_arr_dep(1),kep_arr_dep(2),kep_arr_dep(3),kep_arr_dep(4),kep_arr_dep(5),kep_arr_dep(6),mu_Sun);
[r_arr_fb, ~] = kep2car(kep_arr_fb(1),kep_arr_fb(2),kep_arr_fb(3),kep_arr_fb(4),kep_arr_fb(5),kep_arr_fb(6),mu_Sun);
[r_arr_arr, ~] = kep2car(kep_arr_arr(1),kep_arr_arr(2),kep_arr_arr(3),kep_arr_arr(4),kep_arr_arr(5),kep_arr_arr(6),mu_Sun);

T_dep = 2*pi*sqrt(kep_dep(1)^3/mu_Sun);
T_fb = 2*pi*sqrt(kep_fb_fb(1)^3/mu_Sun);
T_arr = 2*pi*sqrt(kep_arr_arr(1)^3/mu_Sun);

span_dep = (0:24*3600:T_dep);
span_fb = (0:24*3600:T_fb);
span_arr = (0:24*3600:T_arr);

[~,orbit_dep] = ode45(@(t,s) ode_2bp(t,s,mu_Sun),span_dep,[r_dep_Planet; v_dep_Planet],ode_options);
orbit_dep = orbit_dep/AU;
[~,orbit_fb] = ode45(@(t,s) ode_2bp(t,s,mu_Sun),span_fb,[r_fb_dep; v_fb_dep],ode_options);
orbit_fb = orbit_fb/AU;
[~,orbit_arr] = ode45(@(t,s) ode_2bp(t,s,mu_Sun),span_arr,[r_arr_dep; v_arr_dep],ode_options);
orbit_arr = orbit_arr/AU;

ToF_1 = (fbTime-depTime)*24*3600;
span_1 = (0:24*3600:ToF_1);
[~,~,~,~,vdep_t,vfb_t1,~,~] = lambertMR(r_dep_Planet,r_fb_fb,ToF_1,mu_Sun,orbitType,Nrev,Ncase,options);
[~,transfer1] = ode45(@(t,s) ode_2bp(t,s,mu_Sun),span_1,[r_dep_Planet; vdep_t'],ode_options);
transfer1 = transfer1/AU;


ToF_2 = (arrTime-fbTime)*24*3600;
span_2 = (0:24*3600:ToF_2);
[~,~,~,~,vfb_t2,varr_t,~,~] = lambertMR(r_fb_fb,r_arr_arr,ToF_2,mu_Sun,orbitType,Nrev,Ncase,options);
[~,transfer2] = ode45(@(t,s) ode_2bp(t,s,mu_Sun),span_2,[r_fb_fb; vfb_t2'],ode_options);
transfer2 = transfer2/AU;


orbit_dep_rep = repmat(orbit_dep,10,1);
orbit_dep1 = orbit_dep_rep(1:length(span_1),:);
orbit_dep2 = orbit_dep_rep(length(span_1):length(span_1)+length(span_2),:);

orbit_fb_rep = repmat(orbit_fb,6,1);
orbit_fb1 = orbit_fb_rep(1:length(span_1),:);
orbit_fb2 = orbit_fb_rep(length(span_1):length(span_1)+length(span_2),:);

orbit_arr_rep = repmat(orbit_arr,4,1);
orbit_arr1 = orbit_arr_rep(1:length(span_1),:);
orbit_arr2 = orbit_arr_rep(length(span_1):length(span_1)+length(span_2),:);


if nargin == 9
    Nrev = ga2option(1);
    Ncase = ga2option(2);
    ga2Time  = date2mjd2000(ga2Date);
    kep_dep_ga2 = uplanet( ga2Time, dep_Planet);
    kep_fb_ga2  = uplanet( ga2Time, fb_Planet);
    kep_arr_ga2 = uplanet( ga2Time, arr_Planet);
    [r_dep_ga2, ~] = kep2car(kep_dep_ga2(1),kep_dep_ga2(2),kep_dep_ga2(3),kep_dep_ga2(4),kep_dep_ga2(5),kep_dep_ga2(6),mu_Sun);
    [r_fb_ga2, ~] = kep2car(kep_fb_ga2(1),kep_fb_ga2(2),kep_fb_ga2(3),kep_fb_ga2(4),kep_fb_ga2(5),kep_fb_ga2(6),mu_Sun);
    [r_arr_ga2, ~] = kep2car(kep_arr_ga2(1),kep_arr_ga2(2),kep_arr_ga2(3),kep_arr_ga2(4),kep_arr_ga2(5),kep_arr_ga2(6),mu_Sun);
    ToF_ga2 = (depTime-ga2Time)*24*3600;
    span_ga2 = (0:24*3600:ToF_ga2);
    [~,~,~,~,v1_ga2,~,~,~] = lambertMR(r_dep_ga2,r_dep_Planet,ToF_ga2,mu_Sun,orbitType,Nrev,Ncase,options);
    [~,transfer_ga2] = ode45(@(t,s) ode_2bp(t,s,mu_Sun),span_ga2,[r_dep_ga2; v1_ga2'],ode_options);
    transfer_ga2 = transfer_ga2/AU;
    orbit_dep_ga2 = [orbit_dep_rep(end-length(span_ga2)+2:end,:);orbit_dep_rep(1,:)];
    orbit_fb_ga2 = [orbit_fb_rep(end-length(span_ga2)+2:end,:);orbit_fb_rep(1,:)];
    orbit_arr_ga2 = [orbit_arr_rep(end-length(span_ga2)+2:end,:);orbit_arr_rep(1,:)];
else
    %Plot
    Origin=[0 0 0];

    DV_dep_vec=vdep_t-v_dep_planet;
    [p1, p2, p3, p4]=Maneuver_point(v_dep_planet, vdep_t, Origin);
    lgd1=sprintf('Mercury Velocity = %0.5g Km/s', norm(varr_t));
    lgd2=sprintf('V-T1-i = %0.5g Km/s', norm(v_arr_planet));
    lgd3=sprintf('Maneuver Cost = %0.5g Km/s', norm(DV_dep_vec));
    legend( [p1, p2, p3, p4], 'first maneuver point', string(lgd1), string(lgd2), string(lgd3),'Interpreter','latex')
    title('Departure Velocity Triangle')

    DV_arr_vec=v_arr_planet-varr_t;
    [p1, p2, p3, p4]= Maneuver_point(varr_t, v_arr_planet, Origin);
    lgd1=sprintf('V-T2-f = %0.5g Km/s', norm(varr_t));
    lgd2=sprintf('Earth Velocity = %0.5g Km/s', norm(v_arr_planet));
    lgd3=sprintf('Maneuver Cost = %0.5g Km/s', norm(DV_arr_vec));
    legend( [p1, p2, p3, p4], 'final maneuver point', string(lgd1), string(lgd2), string(lgd3))
    title('Arrival Velocity Triangle')

    V_minus=vfb_t1;
    V_plus=vfb_t2;
    DV_fb_vec=V_plus-V_minus;
    [p1, p2, p3, p4]= Maneuver_point(V_minus, V_plus, Origin);
    lgd1=sprintf('V minus = %0.5g Km/s', norm(V_minus));
    lgd2=sprintf('V plus = %0.5g Km/s', norm(V_plus));
    lgd3=sprintf('DV_f_b = %0.5g Km/s', norm(DV_fb_vec));
    legend( [p1, p2, p3, p4], 'GA planet', string(lgd1), string(lgd2), string(lgd3))
    title('Fly-By Velocity Triangle -- Heliocentric ref frame')


    [v_inf_i, v_inf_f, DV_ga_vec, Time_in_SOI, delta_deg, hp]=Hyperbola_plot(V_minus, V_plus, fb_Planet, fbDate);
    [p1, p2, p3, p4, p5]=Maneuver_point(v_inf_i, v_inf_f, Origin, DV_ga_vec);
    lgd1=sprintf('V inf_- = %0.5g Km/s', norm(V_minus));
    lgd2=sprintf('V inf_+ = %0.5g Km/s', norm(V_plus));
    lgd3=sprintf('DV_f_b = %0.5g Km/s', norm(DV_fb_vec));
    lgd4=sprintf('Real Cost: DV_g_a= %0.5g Km/s', norm(DV_ga_vec));
    legend( [p1, p2, p3, p4, p5] , ' maneuver point', string(lgd1), string(lgd2), string(lgd3), string(lgd4))
    title('Fly-By Velocity Triangle and Real Cost -- Venus ref frame')
end


% Plot
switch option
    case 1
        
        figure()
        hold on
        grid on

        p1=plot3(transfer1(:,1),transfer1(:,2),transfer1(:,3),'k', 'LineWidth', 1.2);
        plot3(transfer2(:,1),transfer2(:,2),transfer2(:,3),'k', 'LineWidth', 1.2)
        p2=plot3(orbit_dep(:,1),orbit_dep(:,2),orbit_dep(:,3),'b');
        p3=plot3(orbit_fb(:,1),orbit_fb(:,2),orbit_fb(:,3), 'r');
        p4=plot3(orbit_arr(:,1),orbit_arr(:,2),orbit_arr(:,3),'g');
        axis equal
        scatter3(0,0,0,100,[0.9290 0.6940 0.1250],'filled')

        scatter3(r_dep_Planet(1)/AU,r_dep_Planet(2)/AU,r_dep_Planet(3)/AU,50,'b','filled');
        scatter3(r_dep_fb(1)/AU,r_dep_fb(2)/AU,r_dep_fb(3)/AU,50,'b','filled');
        scatter3(r_dep_arr(1)/AU,r_dep_arr(2)/AU,r_dep_arr(3)/AU,50,'b','filled');


        scatter3(r_fb_dep(1)/AU,r_fb_dep(2)/AU,r_fb_dep(3)/AU,50, 'r','filled')
        scatter3(r_fb_Planet(1)/AU,r_fb_Planet(2)/AU,r_fb_Planet(3)/AU,50, 'r','filled')
        scatter3(r_fb_arr(1)/AU,r_fb_arr(2)/AU,r_fb_arr(3)/AU,50, 'r','filled')

        scatter3(r_arr_dep(1)/AU,r_arr_dep(2)/AU,r_arr_dep(3)/AU,50,'g','filled')
        scatter3(r_arr_fb(1)/AU,r_arr_fb(2)/AU,r_arr_fb(3)/AU,50,'g','filled')
        scatter3(r_arr_Planet(1)/AU,r_arr_Planet(2)/AU,r_arr_Planet(3)/AU,50,'g','filled')

        if nargin == 9
            scatter3(r_dep_ga2(1)/AU,r_dep_ga2(2)/AU,r_dep_ga2(3)/AU,'b','filled')
            plot3(transfer_ga2(:,1),transfer_ga2(:,2),transfer_ga2(:,3),'k')
        end
        xlabel('AU')
        ylabel('AU')
        title('Mission plot')
        legend([p1 p2 p3 p4], 'Interplanetary trajectory','Mercury Orbit','Venus Orbit','Earth Orbit')
        
    case 2
        
        figure()
        hold on
        grid on
                        
        Sc = plot3(nan,nan,nan,'.','Color','k','MarkerSize',30);
        dep = plot3(nan,nan,nan,'.','Color','b','MarkerSize',30);
        fb = plot3(nan,nan,nan,'.','Color','r','MarkerSize',30);
        arr = plot3(nan,nan,nan,'.','Color','g','MarkerSize',30);
        
        plot3(orbit_dep(:,1),orbit_dep(:,2),orbit_dep(:,3),'--b')
        plot3(orbit_fb(:,1),orbit_fb(:,2),orbit_fb(:,3),'--r')
        plot3(orbit_arr(:,1),orbit_arr(:,2),orbit_arr(:,3),'--g')
        axis equal
        scatter3(0,0,0,50,[0.9290 0.6940 0.1250],'filled')
        
        if nargin == 9
        plot3(transfer_ga2(:,1),transfer_ga2(:,2),transfer_ga2(:,3),'k')
        for i = 1:length(span_ga2)
            set(dep,'XData',orbit_dep_ga2(i,1),'YData',orbit_dep_ga2(i,2),'ZData',orbit_dep_ga2(i,3));
            set(fb,'XData',orbit_fb_ga2(i,1),'YData',orbit_fb_ga2(i,2),'ZData',orbit_fb_ga2(i,3));
            set(arr,'XData',orbit_arr_ga2(i,1),'YData',orbit_arr_ga2(i,2),'ZData',orbit_arr_ga2(i,3));
            set(Sc,'XData',transfer_ga2(i,1),'YData',transfer_ga2(i,2),'ZData',transfer_ga2(i,3));
            drawnow
        end
        end
        
        plot3(transfer1(:,1),transfer1(:,2),transfer1(:,3),'k')
        for i = 1:length(span_1)
            set(dep,'XData',orbit_dep1(i,1),'YData',orbit_dep1(i,2),'ZData',orbit_dep1(i,3));
            set(fb,'XData',orbit_fb1(i,1),'YData',orbit_fb1(i,2),'ZData',orbit_fb1(i,3));
            set(arr,'XData',orbit_arr1(i,1),'YData',orbit_arr1(i,2),'ZData',orbit_arr1(i,3));
            set(Sc,'XData',transfer1(i,1),'YData',transfer1(i,2),'ZData',transfer1(i,3));
            drawnow
        end
        
        plot3(transfer2(:,1),transfer2(:,2),transfer2(:,3),'k')
        for i = 1:length(span_2)
            set(dep,'XData',orbit_dep2(i,1),'YData',orbit_dep2(i,2),'ZData',orbit_dep2(i,3));
            set(fb,'XData',orbit_fb2(i,1),'YData',orbit_fb2(i,2),'ZData',orbit_fb2(i,3));
            set(arr,'XData',orbit_arr2(i,1),'YData',orbit_arr2(i,2),'ZData',orbit_arr2(i,3));
            set(Sc,'XData',transfer2(i,1),'YData',transfer2(i,2),'ZData',transfer2(i,3));
            drawnow
        end
        xlabel('AU')
        ylabel('AU')
        title('Mission plot')
end    
        
end
