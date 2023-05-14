%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ASSIGNEMENT 1 - ORBITAL MECHANICS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interplanetary Mission Script

% Mission constraints as indicated in the Assignement:
% Departure Planet:     Earth 
% Fly-by Planet:        Saturn
% Arrival Neo ID:       51
% Earliest Departire:   30-10-26 
% Latest Arrival:       27-04-61

clear
clc

addpath('..\MATLAB functions');              % Add WeBeep Functions - OMIT IN PROJECT DELIVERY
addpath('..\MATLAB functions\time');         % Add WeBeep Functions - OMIT IN PROJECT DELIVERY
addpath('functions\')                        % Add Student Functions
addpath('functions\Data_Images\')            % Add Precomputed Data and Planet Textures
               
%% Time Data

time_departure_earliest = date2mjd2000([2026,10,30,0,0,0]);     % 2026 October  30 00:00:00
time_departure_latest   = date2mjd2000([2061,4,27,0,0,0]);      % 2061 April    27 00:00:00

time_arrival_earliest   = date2mjd2000([2026,10,30,0,0,0]);     % 2026 October  30 00:00:00
time_arrival_latest     = date2mjd2000([2061,4,27,0,0,0]);      % 2061 April    27 00:00:00

%% Celestial Body data

AU          = astroConstants(2);    % Astronimical Unit         [km]

% Planet 1 == 3; EARTH
ID1=3;      % Earth
[Planet_1_kep,ksun]=uplanet(time_departure_earliest,ID1);
mu_Earth=astroConstants(13);
R_earth=astroConstants(23);
Period_Earth=2*pi*sqrt((Planet_1_kep(1)^3)/ksun);

% Planet 2 == 6; SATURN
ID2=6;      % Saturn
[Planet_2_kep,~]=uplanet(time_departure_earliest,ID2);
mu_Saturn=astroConstants(16);
R_saturn=astroConstants(26);
Period_Saturn=2*pi*sqrt((Planet_2_kep(1)^3)/ksun);

% NEO 51 (40 within function)
IDNEO=51;   % NEO 
[NEO_kep,~,~,~]=ephNEO(time_arrival_latest,IDNEO);
Period_NEO=2*pi*sqrt((NEO_kep(1)^3)/ksun);

% Synodic Periods

% Earth to Saturn
Syn_Earth2Saturn=abs(1/((1/Period_Earth)-(1/Period_Saturn)));

% Earth to NEO
Syn_Earth2NEO=abs(1/((1/Period_Earth)-(1/Period_NEO)));

% Saturn to NEO
Syn_Saturn2NEO=abs(1/((1/Period_Saturn)-(1/Period_NEO)));

%% Plot of Orbits - Heliocentric Reference System

time_steps=400;
%time_departure_span     = linspace(date2mjd2000([2026,4,27,0,0,0]),date2mjd2000([2030,1,1,0,0,0]),time_steps);
time_departure_span     = linspace(time_departure_earliest,time_departure_latest,time_steps);

time_arrival_span       = linspace(time_arrival_earliest,time_arrival_latest,time_steps);
%time_arrival_span       = linspace(date2mjd2000([2026,4,27,0,0,0]),date2mjd2000([2030,4,27,0,0,0]),time_steps);


% EARTH
[ys_Planet1_earliest_departure,v_Planet1_earliest_departure]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

s_Planet1_earliest_departure = [ys_Planet1_earliest_departure';v_Planet1_earliest_departure'];
T_Planet1 = 2*pi*sqrt( Planet_1_kep(1)^3/ksun );
t_Planet1_orbit_t_span=linspace(time_departure_earliest*24*3600,time_departure_earliest*24*3600+T_Planet1,5000);

[~,ys_Planet1_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet1_orbit_t_span,s_Planet1_earliest_departure);

% SATURN
[ys_Planet2_earliest_departure,v_Planet2_earliest_departure]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

s_Planet2_earliest_departure = [ys_Planet2_earliest_departure';v_Planet2_earliest_departure'];
T_Planet2 = 2*pi*sqrt( Planet_2_kep(1)^3/ksun );
t_Planet2_orbit_t_span=linspace(time_departure_earliest*24*3600,time_departure_earliest*24*3600+T_Planet2,5000);

[~,ys_Planet2_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet2_orbit_t_span,s_Planet2_earliest_departure);

% NEO
[ys_NEO_earliest_departure,v_NEO_earliest_departure]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

s_NEO_earliest_departure = [ys_NEO_earliest_departure';v_NEO_earliest_departure'];
T_NEO = 2*pi*sqrt( NEO_kep(1)^3/ksun );
t_NEO_orbit_t_span=linspace(time_departure_earliest*24*3600,time_departure_earliest*24*3600+T_NEO,5000);

[~,ys_NEO_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_NEO_orbit_t_span,s_NEO_earliest_departure);

% Plotting - HELIOCENTRIC - AU
figure(1)

% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/3;
Yt=Yt/3;
Zt=Zt/3;
image_file = 'Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

title('Orbit of Earth, Saturn and NEO','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color','b')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color','r')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color','g')

legend('Earth Orbit',...
    'Saturn Orbit',...
    'NEO Orbit','Interpreter','latex','Fontsize',11)

%% Lambert Solution for Earth to Saturn - COMPLETE TIME SPAN

deltaV=zeros(time_steps);
h = waitbar(0,'Please wait...');
tic;
for i=1:time_steps % departure
    
    for j=1:time_steps % arrival

        if time_arrival_span(j)<time_departure_span(i)

            deltaV(i,j)=NaN;
        else
            [deltaV_it,~,~,~,~]=delta_V_Interplanetary(time_departure_span(i),time_arrival_span(j),ID1,ID2);
            
            deltaV(i,j)=deltaV_it;
        end
    end
    
    waitbar(i / time_steps)

end
toc;

%% CONTOUR PLOT - Earth to Saturn - COMPLETE TIME SPAN
[time_departure_mesh,time_arrival_mesh]=meshgrid(time_departure_span,time_arrival_span);
figure(2)
% contour(time_departure_mesh,time_arrival_mesh,deltaV',[15,20,25,30],'ShowText','on');
contour(time_departure_mesh,time_arrival_mesh,deltaV',[15,20,25,30,40,80,100]);
colorbar
% xticks(linspace(time_departure_earliest,time_departure_latest,5))
% xticklabels({"01/04/2003","01/05/2003","01/06/2003","01/07/2003","01/08/2003"})
% yticks(linspace(time_arrival_earliest,time_arrival_latest,7))
% yticklabels({"01/09/2003","01/10/2003","01/11/2003","01/12/2003","01/01/2004","01/02/2004","01/03/2004"})
xlabel('Departure Time')
ylabel('Arrival Time')
grid on

%% Lambert Solution for Saturn to NEO - COMPLETE TIME SPAN

deltaV2=zeros(time_steps);
h = waitbar(0,'Please wait...');

for i=1:time_steps % departure
    
    for j=1:time_steps % arrival

        if time_arrival_span(j)<time_departure_span(i)

            deltaV2(i,j)=NaN;
        else
            [deltaV_it,~,~,~,~]=delta_V_Interplanetary(time_departure_span(i),time_arrival_span(j),ID2,IDNEO);

            deltaV2(i,j)=deltaV_it;
        end

    end
    
    waitbar(i / time_steps)

end

%% CONTOUR PLOT - Saturn to NEO - COMPLETE TIME SPAN
[time_departure_mesh,time_arrival_mesh]=meshgrid(time_departure_span,time_arrival_span);
figure(3)
% contour(time_departure_mesh,time_arrival_mesh,deltaV',[15,20,25,30],'ShowText','on');
contour(time_departure_mesh,time_arrival_mesh,deltaV2',[10,15,20,25,30,40,60,80]);
colorbar
% xticks(linspace(time_departure_earliest,time_departure_latest,5))
% xticklabels({"01/04/2003","01/05/2003","01/06/2003","01/07/2003","01/08/2003"})
% yticks(linspace(time_arrival_earliest,time_arrival_latest,7))
% yticklabels({"01/09/2003","01/10/2003","01/11/2003","01/12/2003","01/01/2004","01/02/2004","01/03/2004"})
xlabel('Departure Time')
ylabel('Arrival Time')
grid on

%% COARSE COMPUTATION
ToF_GA_vector       = 300:30:12450; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
ToF_arrival_vector  = 100:30:12650; % DAYS - MUST CONVERT IN SECONDS for LAMBERT

len_GA      =length(ToF_GA_vector);
len_arrival =length(ToF_arrival_vector);

% MATRICES
deltaV_matrix=NaN(time_steps,len_GA,len_arrival);

deltaV_GA_matrix=NaN(time_steps,len_GA,len_arrival);
deltaV_departure_matrix=NaN(time_steps,len_GA,len_arrival);
deltaV_arrival_matrix=NaN(time_steps,len_GA,len_arrival);

rp_matrix=NaN(time_steps,len_GA,len_arrival);  

fsolve_error_matrix=NaN(time_steps,len_GA,len_arrival);

delta_angle_matrix=NaN(time_steps,len_GA,len_arrival);

rp_crit=R_saturn*1.3; % Guess. 30% extra from Saturn Radius

h = waitbar(0,'Please wait...');

tic;
for i=1:time_steps % Departure Date

    time_departure=time_departure_span(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Position and Velocity of Earth
    [Earth_kep,~]=uplanet(time_departure,ID1);

    [ys_Earth_tdeparture,~]=kep2car(Earth_kep(1),...
        Earth_kep(2),...
        rad2deg(Earth_kep(3)),...
        rad2deg(Earth_kep(4)),...
        rad2deg(Earth_kep(5)),...
        rad2deg(Earth_kep(6)),...
        ksun);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for j=1:len_GA % ToF Gravity Assist
        
        ToF_GA  = ToF_GA_vector(j);
        time_GA = time_departure+ToF_GA;
        
        if time_GA > time_arrival_latest
            break
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Position and Velocity of Saturn
        [Saturn_kep,~]=uplanet(time_GA,ID2);

        [ys_Saturn_tGA,V_Saturn_tGA]=kep2car(Saturn_kep(1),...
            Saturn_kep(2),...
            rad2deg(Saturn_kep(3)),...
            rad2deg(Saturn_kep(4)),...
            rad2deg(Saturn_kep(5)),...
            rad2deg(Saturn_kep(6)),...
            ksun);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ToF convert from days to seconds

        % Final Velocity from Lambert is Velocity of Arrival at Saturn
        [~,deltaV_departure,~,~,V_inf_minus]=delta_V_Interplanetary(time_departure,time_GA,ID1,ID2);

        % Cutoff deltaV manoeuvre from Earth to transfer arc to Saturn
        if norm(deltaV_departure)>16
            continue
        end

        for k=1:len_arrival

            ToF_arrival     = ToF_arrival_vector(k); % DAYS
            time_arrival    = time_GA+ToF_arrival; % DAYS
                
            if time_arrival > time_arrival_latest
                break
            end

            %error_matrix(i,j,k)=Error_lambert;

            [NEO_kep, ~, ~] = ephNEO(time_arrival, IDNEO);

            [ys_NEO_tarrival,~]=kep2car(NEO_kep(1),...
                NEO_kep(2),...
                rad2deg(NEO_kep(3)),...
                rad2deg(NEO_kep(4)),...
                rad2deg(NEO_kep(5)),...
                rad2deg(NEO_kep(6)),...
                ksun);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Obtain V_inf_plus to do the transfer arc. Equivalent to
            % Initial Velocity in Lambert Solution
        
            % Convert ToF from DAYS to SECONDS
            [~,~,V_inf_plus,deltaV_arrival,~]=delta_V_Interplanetary(time_GA,time_arrival,ID2,IDNEO);
            %[~,~,~,Error_lambert,V_inf_plus,~,~,~] = lambertMR( ys_Saturn_tGA, ys_NEO_tarrival, ToF_arrival*(24*3600), ksun, 0, 0, 0 );
            
            if norm(deltaV_arrival)>10
                continue
            end

            v_inf_minus = V_inf_minus-V_Saturn_tGA;     % Saturn's Point of view
            v_inf_plus  = V_inf_plus-V_Saturn_tGA;

            delta_angle = acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_minus)*norm(v_inf_plus)));     % Total turn angle is the angle between v_inf_minus and v_inf_plus

            delta_angle_matrix(i,j,k)=delta_angle;

            rp_specific = @(rp) rp_solver(rp,norm(v_inf_minus),norm(v_inf_plus),delta_angle,mu_Saturn);

            [rp,~,rp_sol_flag,~] = fzero(rp_specific,R_saturn*10,optimset('Display','off'));
            
            %rp=real(rp);
            fsolve_error_matrix(i,j,k)=rp_sol_flag;
            rp_matrix(i,j,k)=rp;

            if (rp<rp_crit) || (rp_sol_flag<1) % error flag when 0,-1,-2,-3,-4
                continue
            end            

            v_inf_minus_sq_norm = sum(v_inf_minus.*v_inf_minus); % Norm of vector squared. Reduces computational cost
            v_inf_plus_sq_norm  = sum(v_inf_plus.*v_inf_plus);

            vp_minus_norm     = sqrt(((rp*v_inf_minus_sq_norm+2*mu_Saturn)/rp));  % Perifocal Velocity Magnitude at Perigee from Incoming Trajectory [km/s]
            vp_plus_norm      = sqrt(((rp*v_inf_plus_sq_norm+2*mu_Saturn)/rp));   % Perifocal Velocity Magnitude at Perigee from Outgoing Trajectory [km/s]

            % Powered fly-by deltaV manoeuvre
            deltaV_GA   = abs(vp_plus_norm-vp_minus_norm);

            deltaV_matrix(i,j,k)=deltaV_GA+norm(deltaV_departure)+norm(deltaV_arrival);
            % debugging
            deltaV_GA_matrix(i,j,k)=deltaV_GA;
            deltaV_departure_matrix(i,j,k)=norm(deltaV_departure);
            deltaV_arrival_matrix(i,j,k)=norm(deltaV_arrival);


        end
    end
    waitbar(i/time_steps)
end

num_computed=sum(sum(sum(~isnan(deltaV_matrix))))*100/(time_steps*len_arrival*len_GA); % Percentage of computed deltaV Values

[minDeltaV, idx] = min(deltaV_matrix(:));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

% min DeltaV
fprintf('Minimum deltaV: %g \nDeltaV departure: %g \nDeltaV GA: %g \nDeltaV arrival: %g \n',minDeltaV,deltaV_departure_matrix(n,m,t),deltaV_GA_matrix(n,m,t),deltaV_arrival_matrix(n,m,t))

time_departure_minimum=mjd20002date(time_departure_span(n));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m)+ToF_arrival_vector(t));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************************** \n');
fprintf('Radius of periguee: %g x Radius of Saturn\nTurn Angle (deg): %g\nToF to Saturn (days): %g\nToF to NEO (days): %g\n',rp_matrix(n,m,t)/R_saturn,rad2deg(delta_angle_matrix(n,m,t)),ToF_GA_vector(m),ToF_arrival_vector(t))
fprintf('Computed Percentage: %g\n',num_computed);

%% Plot Heliocentric Orbit from best coarse solution

load('first_coarse.mat')
[~, idx] = min(deltaV_matrix(:));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

% Times - Coarse
best_coarse_departure=time_departure_span(n); %  mjd2000
best_coarse_time_GA=time_departure_span(n)+ToF_GA_vector(m); 
best_coarse_time_arrival=best_coarse_time_GA+ToF_arrival_vector(t);

% EARTH
[Planet_1_kep,ksun]=uplanet(best_coarse_departure,ID1);
[ys_Planet1_coarse_departure,v_Planet1_coarse_departure]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

s_Planet1_coarse_departure = [ys_Planet1_coarse_departure';v_Planet1_coarse_departure'];
T_Planet1 = 2*pi*sqrt( Planet_1_kep(1)^3/ksun );
t_Planet1_orbit_t_span=linspace(best_coarse_departure*24*3600,best_coarse_departure*24*3600+T_Planet1,5000);

[~,ys_Planet1_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet1_orbit_t_span,s_Planet1_coarse_departure);

% SATURN
[Planet_2_kep,ksun]=uplanet(best_coarse_time_GA,ID2);
[ys_Planet2_coarse_GA,v_Planet2_coarse_GA]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

s_Planet2_coarse_GA = [ys_Planet2_coarse_GA';v_Planet2_coarse_GA'];
T_Planet2 = 2*pi*sqrt( Planet_2_kep(1)^3/ksun );
t_Planet2_orbit_t_span=linspace(best_coarse_time_GA*24*3600,best_coarse_time_GA*24*3600+T_Planet2,5000);

[~,ys_Planet2_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet2_orbit_t_span,s_Planet2_coarse_GA);

% NEO
[NEO_kep,~,~,~]=ephNEO(best_coarse_time_arrival,IDNEO);
[ys_NEO_coarse_arrival,v_NEO_coarse_arrival]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

s_NEO_coarse_arrival = [ys_NEO_coarse_arrival';v_NEO_coarse_arrival'];
T_NEO = 2*pi*sqrt( NEO_kep(1)^3/ksun );
t_NEO_orbit_t_span=linspace(best_coarse_time_arrival*24*3600,best_coarse_time_arrival*24*3600+T_NEO,5000);

[~,ys_NEO_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_NEO_orbit_t_span,s_NEO_coarse_arrival);

% Transfer Arc - Earth to Saturn
ToF_GA_span=linspace(best_coarse_departure*24*3600,(best_coarse_departure+ToF_GA_vector(m))*24*3600,1000);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet1_coarse_departure, ys_Planet2_coarse_GA, ToF_GA_vector(m)*24*3600, ksun, 0, 0, 0 );
s_transfer_GA=[ys_Planet1_coarse_departure',VI'];
[~,ys_trans_GA] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_GA_span,s_transfer_GA);


% Transfer Arc - Saturn to NEO
ToF_arrival_span=linspace(best_coarse_time_GA*24*3600,(best_coarse_time_GA+ToF_arrival_vector(t))*24*3600,1000);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet2_coarse_GA, ys_NEO_coarse_arrival, ToF_arrival_vector(t)*24*3600, ksun, 0, 0, 0 );
s_transfer_arrival=[ys_Planet2_coarse_GA',VI'];
[~,ys_trans_arrival] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_arrival_span,s_transfer_arrival);

% Plotting - HELIOCENTRIC - AU
figure(50)

% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/3;
Yt=Yt/3;
Zt=Zt/3;
image_file = 'https://upload.wikimedia.org/wikipedia/commons/9/99/Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

title('Orbit of Earth, Saturn and NEO','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color','b')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color','r')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color','g')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color','m','LineStyle','--')

% Earth to Saturn arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color','c','LineStyle','--')

legend('Earth Orbit',...
    'Saturn Orbit',...
    'NEO Orbit',...
    'Saturn Transfer',...
    'NEO transfer','Interpreter','latex','Fontsize',11)


%% Plot - Specific deltaV Index Selection

deltaV_array=sort(deltaV_matrix(deltaV_matrix>0));

% Start in N best deltaV
idx2=find(deltaV_matrix==deltaV_array(1));
[n2, m2, t2] = ind2sub(size(deltaV_matrix),idx2);

% Times - Coarse
best_coarse_departure=time_departure_span(n2); %  mjd2000
best_coarse_time_GA=time_departure_span(n2)+ToF_GA_vector(m2); 
best_coarse_time_arrival=best_coarse_time_GA+ToF_arrival_vector(t2);

% EARTH
[Planet_1_kep,ksun]=uplanet(best_coarse_departure,ID1);
[ys_Planet1_coarse_departure,v_Planet1_coarse_departure]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

s_Planet1_coarse_departure = [ys_Planet1_coarse_departure';v_Planet1_coarse_departure'];
T_Planet1 = 2*pi*sqrt( Planet_1_kep(1)^3/ksun );
t_Planet1_orbit_t_span=linspace(best_coarse_departure*24*3600,best_coarse_departure*24*3600+T_Planet1,5000);

[~,ys_Planet1_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet1_orbit_t_span,s_Planet1_coarse_departure);

% SATURN
[Planet_2_kep,ksun]=uplanet(best_coarse_time_GA,ID2);
[ys_Planet2_coarse_GA,v_Planet2_coarse_GA]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

s_Planet2_coarse_GA = [ys_Planet2_coarse_GA';v_Planet2_coarse_GA'];
T_Planet2 = 2*pi*sqrt( Planet_2_kep(1)^3/ksun );
t_Planet2_orbit_t_span=linspace(best_coarse_time_GA*24*3600,best_coarse_time_GA*24*3600+T_Planet2,5000);

[~,ys_Planet2_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet2_orbit_t_span,s_Planet2_coarse_GA);

% NEO
[NEO_kep,~,~,~]=ephNEO(best_coarse_time_arrival,IDNEO);
[ys_NEO_coarse_arrival,v_NEO_coarse_arrival]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

s_NEO_coarse_arrival = [ys_NEO_coarse_arrival';v_NEO_coarse_arrival'];
T_NEO = 2*pi*sqrt( NEO_kep(1)^3/ksun );
t_NEO_orbit_t_span=linspace(best_coarse_time_arrival*24*3600,best_coarse_time_arrival*24*3600+T_NEO,5000);

[~,ys_NEO_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_NEO_orbit_t_span,s_NEO_coarse_arrival);

% Transfer Arc - Earth to Saturn
ToF_GA_span=linspace(best_coarse_departure*24*3600,(best_coarse_departure+ToF_GA_vector(m2))*24*3600,1000);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet1_coarse_departure, ys_Planet2_coarse_GA, ToF_GA_vector(m2)*24*3600, ksun, 0, 0, 0 );
s_transfer_GA=[ys_Planet1_coarse_departure',VI'];
[~,ys_trans_GA] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_GA_span,s_transfer_GA);


% Transfer Arc - Saturn to NEO
ToF_arrival_span=linspace(best_coarse_time_GA*24*3600,(best_coarse_time_GA+ToF_arrival_vector(t2))*24*3600,1000);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet2_coarse_GA, ys_NEO_coarse_arrival, ToF_arrival_vector(t2)*24*3600, ksun, 0, 0, 0 );
s_transfer_arrival=[ys_Planet2_coarse_GA',VI'];
[~,ys_trans_arrival] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_arrival_span,s_transfer_arrival);


% Plotting - HELIOCENTRIC - AU
figure(51)

% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/3;
Yt=Yt/3;
Zt=Zt/3;
image_file = 'https://upload.wikimedia.org/wikipedia/commons/9/99/Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

title('Orbit of Earth, Saturn and NEO','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color','b')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color','r')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color','g')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color','m','LineStyle','-')

% Earth to Saturn arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color','c','LineStyle','-')

legend('Earth Orbit',...
    'Saturn Orbit',...
    'NEO Orbit',...
    'Saturn Transfer',...
    'NEO transfer','Interpreter','latex','Fontsize',11)

fprintf('N-Minimum deltaV: %g \nDeltaV departure: %g \nDeltaV GA: %g \nDeltaV arrival: %g \n',deltaV_matrix(n2,m2,t2),deltaV_departure_matrix(n2,m2,t2),deltaV_GA_matrix(n2,m2,t2),deltaV_arrival_matrix(n2,m2,t2))

time_departure_minimum=mjd20002date(time_departure_span(n2));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(time_departure_span(n2)+ToF_GA_vector(m2));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(time_departure_span(n2)+ToF_GA_vector(m2)+ToF_arrival_vector(t2));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************************** \n');
fprintf('Radius of periguee: %g x Radius of Saturn\nTurn Angle (deg): %g\nToF to Saturn (days): %g\nToF to NEO (days): %g\n',rp_matrix(n2,m2,t2)/R_saturn,rad2deg(delta_angle_matrix(n2,m2,t2)),ToF_GA_vector(m2),ToF_arrival_vector(t2))


%% Post Processing COARSE

load("first_coarse.mat")
time_steps=size(deltaV_matrix,1);

num_computed=sum(sum(sum(~isnan(deltaV_matrix))))*100/(time_steps*size(deltaV_matrix,2)*size(deltaV_matrix,3)); % Percentage of computed deltaV Values


best_limit=20000;

best_deltaV=sort(deltaV_matrix(:)); % Best deltaV
best_deltaV=best_deltaV(1:best_limit);

[~,idx_array]=sort(deltaV_matrix(:)); % CREATES A SORTED COLUMN VECTOR OF EXISTING DELTAV values

[n2, m2, t2] = ind2sub(size(deltaV_matrix),idx_array);

n2=n2(1:best_limit);
m2=m2(1:best_limit);
t2=t2(1:best_limit);

cmap=colormap(jet(best_limit));

hfig=figure(100);
scatter3(time_departure_span(n2),ToF_GA_vector(m2),ToF_arrival_vector(t2),5,cmap,'filled','MarkerFaceAlpha',.4)
clim([min(best_deltaV),max(best_deltaV)])
colorbar

title('Mission deltaV - Coarse - complete time span','Interpreter','latex','Fontsize',11);
xlabel('departure time [days - mjd2000]','Interpreter','latex','Fontsize',11); 
ylabel('Time of Flight to Saturn [days]','Interpreter','latex','Fontsize',11); 
zlabel('Time of Flight to NEO [days]','Interpreter','latex','Fontsize',11);
xlim([time_departure_span(1),time_departure_span(end)])

set(gca,'YLim',[0 8000],'YTick',0:2000:8000)
set(gca,'ZLim',[0 10000],'ZTick',0:2000:10000)
% axis equal;
grid on;
view([-21,28])

fname = 'Coarse_Grid';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

%set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')




% Plot of Time windows with percentage of cases 
iscomputed=~isnan(deltaV_matrix); % Logical 1 is deltaV(i,j,k) exists

departure_computed=zeros(time_steps,2);
for k=1:time_steps
    departure_computed(k,1)=sum(iscomputed(k,:,:),'all');
    if sum(iscomputed(k,:,:),'all')>0
        departure_computed(k,2)=1;
    end
end

hfig=figure(101);
set(hfig,'defaultAxesColorOrder',[[0, 0, 0]; [0, 0, 0]]);
yyaxis left
plot(time_departure_span,departure_computed(:,1),'LineWidth',1.5)
title('Best $\Delta v$ - Coarse - complete time span','Interpreter','latex','Fontsize',11);
xlabel('departure time [days - mjd2000]','Interpreter','latex','Fontsize',11); 
ylabel('Number of Cases','Interpreter','latex','Fontsize',11); 
yyaxis right
plot(time_departure_span,departure_computed(:,1)*100/sum(iscomputed,'all'),'LineStyle','none','Marker','none')
ylabel('Percentage of Total Cases [\%]','Interpreter','latex','Fontsize',11); 

%hfig = figure;  % save the figure handle in a variable
fname = 'Coarse_Number_Cases';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')


%% COARSE COMPUTATION SECOND VERSION - Refined Time window

%long
% ToF_GA_vector       = 2600:10:3800; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% ToF_arrival_vector  = 2400:10:4200; % DAYS - MUST CONVERT IN SECONDS for LAMBERT

% short
ToF_GA_vector       = 900:10:1900; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
ToF_arrival_vector  = 900:10:2200; % DAYS - MUST CONVERT IN SECONDS for LAMBERT

% ToF_GA_vector = 2700:5:2950;
% ToF_arrival_vector = 3300:5:3600;

% ToF_GA_vector=2760;
% ToF_arrival_vector=3450;

time_steps = 200;
%time_departure_span     = linspace(9956.37,10240.5,time_steps);    %1 lon
%time_departure_span     = linspace(10335.3,10578.8,time_steps);    %2  lon 
%time_departure_span     = linspace(10714.1,10966.7,time_steps);    %3 sho+lon
%time_departure_span     = linspace(11093,11345.6,time_steps);     %4  sho
time_departure_span     = linspace(11471.9,11692.9,time_steps);   %5  sho
%time_departure_span     = linspace(11850.8,12071.8,time_steps);   %6  sho
%time_departure_span     = linspace(12229.7,12450.7,time_steps);   %7
%time_departure_span     = linspace(12640.2,12798,time_steps);     %8
%time_departure_span     = linspace(13019,13176.9,time_steps);     %9
%time_departure_span     = linspace(13397.9,13555.8,time_steps);   %10
%time_departure_span     = linspace(13776.8,13903.1,time_steps);   %11
%time_departure_span     = linspace(14155.7,14282,time_steps);     %12
%time_departure_span     = linspace(14534.6,14660.9,time_steps);   %13
%time_departure_span     = linspace(14913.5,15008.2,time_steps);   %14


% Best dVs (less than 18 km/s) where given in departure windows: 1-6

% After second refinement (constrained ToF), best dV: 5, 1, 3, 6, 2, 4
% 5th: 16.9032
% 1st: 17.3419 
% 3rd: 17.4788
% 6th: 17.5593 
% 2nd: 17.5812 
% 4th: 17.6101 

% Third refinement on best window candidates: 5, 1, 3, 6
% time_steps = 20;

% 5th
% ToF_GA_vector       = 1380:2:1460; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% ToF_arrival_vector  = 1530:2:1610; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% time_departure_span     = linspace(11544.1,11557.4,time_steps);   %5  sho

% 1st
% ToF_GA_vector       = 3180:2:3450; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% ToF_arrival_vector  = 2620:2:3270; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% time_departure_span     = linspace(10063.5,10076.3,time_steps);   %5  sho

% 3rd
% ToF_GA_vector       = 1420:2:1670; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% ToF_arrival_vector  = 1360:2:1610; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% time_departure_span     = linspace(10782.6,10804.2,time_steps);   %5  sho

% 6th
% ToF_GA_vector       = 1120:2:1390; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% ToF_arrival_vector  = 1220:2:1490; % DAYS - MUST CONVERT IN SECONDS for LAMBERT
% time_departure_span     = linspace(11916.3,11944.1,time_steps);   %5  sho

len_GA      =length(ToF_GA_vector);
len_arrival =length(ToF_arrival_vector);

% MATRICES
deltaV_matrix=NaN(time_steps,len_GA,len_arrival);

deltaV_GA_matrix=NaN(time_steps,len_GA,len_arrival);
deltaV_departure_matrix=NaN(time_steps,len_GA,len_arrival);
deltaV_arrival_matrix=NaN(time_steps,len_GA,len_arrival);

rp_matrix=NaN(time_steps,time_steps,time_steps); 
r_helio_matrix=NaN(time_steps,time_steps,time_steps); 

error_matrix=NaN(time_steps,len_GA,len_arrival);
fsolve_error_matrix=NaN(time_steps,len_GA,len_arrival);

delta_angle_matrix=NaN(time_steps,len_GA,len_arrival);

rp_crit=R_saturn*1.3; % Guess. 30% extra from Saturn Radius

h = waitbar(0,'Please wait...');

tic;
for i=1:time_steps % Departure Date

    time_departure=time_departure_span(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Position and Velocity of Earth
    [Earth_kep,~]=uplanet(time_departure,ID1);

    [ys_Earth_tdeparture,~]=kep2car(Earth_kep(1),...
        Earth_kep(2),...
        rad2deg(Earth_kep(3)),...
        rad2deg(Earth_kep(4)),...
        rad2deg(Earth_kep(5)),...
        rad2deg(Earth_kep(6)),...
        ksun);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for j=1:len_GA % ToF Gravity Assist
        
        ToF_GA  = ToF_GA_vector(j);
        time_GA = time_departure+ToF_GA;
        
        if time_GA > time_arrival_latest
            break
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Position and Velocity of Saturn
        [Saturn_kep,~]=uplanet(time_GA,ID2);

        [ys_Saturn_tGA,V_Saturn_tGA]=kep2car(Saturn_kep(1),...
            Saturn_kep(2),...
            rad2deg(Saturn_kep(3)),...
            rad2deg(Saturn_kep(4)),...
            rad2deg(Saturn_kep(5)),...
            rad2deg(Saturn_kep(6)),...
            ksun);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ToF convert from days to seconds

        % Final Velocity from Lambert is Velocity of Arrival at Saturn
        [~,deltaV_departure,~,~,V_inf_minus]=delta_V_Interplanetary(time_departure,time_GA,ID1,ID2);

        % Cutoff deltaV manoeuvre from Earth to transfer arc to Saturn
        if norm(deltaV_departure)>13
            continue
        end

        for k=1:len_arrival

            ToF_arrival     = ToF_arrival_vector(k); % DAYS
            time_arrival    = time_GA+ToF_arrival; % DAYS
                
            if time_arrival > time_arrival_latest
                break
            end

            %error_matrix(i,j,k)=Error_lambert;

            [NEO_kep, ~, ~] = ephNEO(time_arrival, IDNEO);

            [ys_NEO_tarrival,~]=kep2car(NEO_kep(1),...
                NEO_kep(2),...
                rad2deg(NEO_kep(3)),...
                rad2deg(NEO_kep(4)),...
                rad2deg(NEO_kep(5)),...
                rad2deg(NEO_kep(6)),...
                ksun);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Obtain V_inf_plus to do the transfer arc. Equivalent to
            % Initial Velocity in Lambert Solution
        
            % Convert ToF from DAYS to SECONDS
            [~,~,V_inf_plus,deltaV_arrival,~]=delta_V_Interplanetary(time_GA,time_arrival,ID2,IDNEO);
            
            if norm(deltaV_arrival)>10
                continue
            end

            v_inf_minus = V_inf_minus-V_Saturn_tGA;     % Saturn's Point of view
            v_inf_plus  = V_inf_plus-V_Saturn_tGA;

            delta_angle = acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_minus)*norm(v_inf_plus)));     % Total turn angle is the angle between v_inf_minus and v_inf_plus

            delta_angle_matrix(i,j,k)=delta_angle;

            rp_specific = @(rp) rp_solver(rp,norm(v_inf_minus),norm(v_inf_plus),delta_angle,mu_Saturn);

            [rp,~,rp_sol_flag,~] = fzero(rp_specific,R_saturn*10,optimset('Display','off'));
            
            %rp=real(rp);
            fsolve_error_matrix(i,j,k)=rp_sol_flag;
            rp_matrix(i,j,k)=rp;

            if (rp<rp_crit) || (rp_sol_flag<1) % error flag when 0,-1,-2,-3,-4
                continue
            end            

            v_inf_minus_sq_norm = sum(v_inf_minus.*v_inf_minus); % Norm of vector squared. Reduces computational cost
            v_inf_plus_sq_norm  = sum(v_inf_plus.*v_inf_plus);

            vp_minus_norm     = sqrt(((rp*v_inf_minus_sq_norm+2*mu_Saturn)/rp));  % Perifocal Velocity Magnitude at Perigee from Incoming Trajectory [km/s]
            vp_plus_norm      = sqrt(((rp*v_inf_plus_sq_norm+2*mu_Saturn)/rp));   % Perifocal Velocity Magnitude at Perigee from Outgoing Trajectory [km/s]

            % Powered fly-by deltaV manoeuvre
            deltaV_GA   = abs(vp_plus_norm-vp_minus_norm);

            deltaV_matrix(i,j,k)=deltaV_GA+norm(deltaV_departure)+norm(deltaV_arrival);
            % debugging
            deltaV_GA_matrix(i,j,k)=deltaV_GA;
            deltaV_departure_matrix(i,j,k)=norm(deltaV_departure);
            deltaV_arrival_matrix(i,j,k)=norm(deltaV_arrival);


        end
    end
    waitbar(i/time_steps)
end

num_computed=sum(sum(sum(~isnan(deltaV_matrix))))*100/(time_steps*len_arrival*len_GA); % Percentage of computed deltaV Values

[minDeltaV, idx] = min(deltaV_matrix(:));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

% min DeltaV
fprintf('Minimum deltaV: %g \nDeltaV departure: %g \nDeltaV GA: %g \nDeltaV arrival: %g \n',minDeltaV,deltaV_departure_matrix(n,m,t),deltaV_GA_matrix(n,m,t),deltaV_arrival_matrix(n,m,t))

time_departure_minimum=mjd20002date(time_departure_span(n));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m)+ToF_arrival_vector(t));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************************** \n');
fprintf('Radius of periguee: %g x Radius of Saturn\nTurn Angle (deg): %g\nToF to Saturn (days): %g\nToF to NEO (days): %g\n',rp_matrix(n,m,t)/R_saturn,rad2deg(delta_angle_matrix(n,m,t)),ToF_GA_vector(m),ToF_arrival_vector(t))


%% Post Processing COARSE SECOND Refinement

% Remove deltaV too high
deltaV_matrix(deltaV_matrix>22)=NaN;

best_limit=1000;

best_deltaV=sort(deltaV_matrix(:)); % Best deltaV
best_deltaV=best_deltaV(1:best_limit);

[~,idx_array]=sort(deltaV_matrix(:)); % CREATES A SORTED COLUMN VECTOR OF EXISTING DELTAV values

[n2, m2, t2] = ind2sub(size(deltaV_matrix),idx_array);

n2=n2(1:best_limit);
m2=m2(1:best_limit);
t2=t2(1:best_limit);

cmap=colormap(jet(best_limit));

hfig=figure(302);
scatter3(time_departure_span(n2),ToF_GA_vector(m2),ToF_arrival_vector(t2),5,cmap,'filled','MarkerFaceAlpha',.7)
clim([min(best_deltaV),max(best_deltaV)])
colorbar

title('5th Departure Window - Mission $\Delta v$ - Refined','Interpreter','latex','Fontsize',11);
%title('Nth Departure Window - Mission $\Delta v$ - Refined','Interpreter','latex','Fontsize',11);
xlabel('departure time [days - mjd2000]','Interpreter','latex','Fontsize',11); 
ylabel('Time of Flight to Saturn','Interpreter','latex','Fontsize',11); 
zlabel('Time of Flight to NEO','Interpreter','latex','Fontsize',11);
%xlim([time_departure_span(1),time_departure_span(end)])
%xlim([time_departure_earliest,time_departure_latest])
% axis equal;
grid on;

view([-44,12])

% fname = 'Refined_Grid_5_Dep';
% 
% picturewidth = 20; % set this parameter and keep it forever
% hw_ratio = 0.65; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% 
% %set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hfig,fname,'-dpdf','-painters','-fillpage')
% %print(hfig,fname,'-dpng','-painters')


% Comparing ToF
hfig = figure(300)
plot(best_deltaV,ToF_GA_vector(m2),'Color',[0.9290 0.6940 0.1250 0.3])
hold on
plot(best_deltaV,ones(1,length(best_deltaV))*mean(ToF_GA_vector(m2)),'Color',[0.9290 0.6940 0.1250 1],'LineWidth',2)
plot(best_deltaV,ToF_arrival_vector(t2),'Color',[0 0.4470 0.7410 0.3])
plot(best_deltaV,ones(1,length(best_deltaV))*mean(ToF_arrival_vector(t2)),'Color',[0 0.4470 0.7410 1],'LineWidth',2)
title('5th Departure Window - ToF comparison','Interpreter','latex','Fontsize',11);
xlabel('Best $\Delta v$ results [km/s]','Interpreter','latex','Fontsize',11); 
ylabel('Time of Flight [days]','Interpreter','latex','Fontsize',11); 
legend('',...
    'ToF to Saturn',...
    '',...
    'ToF to NEO',...
    'Interpreter','latex','Fontsize',11)

fname = 'ToF_comparison_5';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

%set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')


figure(301)
plot(best_deltaV,time_departure_span(n2))

fprintf('In the best %g cases:\nToF_GA range: [%g,%g] days\nToF_arrival range: [%g,%g] days\n',...
    best_limit,min(ToF_GA_vector(m2)),max(ToF_GA_vector(m2)),min(ToF_arrival_vector(t2)),max(ToF_arrival_vector(t2)));

fprintf('Time of Departure: [%g,%g]\n',min(time_departure_span(n2)),max(time_departure_span(n2)))

%% Plot Heliocentric Orbit from best refined solution - Chose iteration 

deltaV_array=sort(deltaV_matrix(deltaV_matrix>0));

% start in N best deltaV
idx2=find(deltaV_matrix==deltaV_array(1));
[n, m, t] = ind2sub(size(deltaV_matrix),idx2);


% Times - Coarse
best_coarse_departure=time_departure_span(n); %  mjd2000
best_coarse_time_GA=time_departure_span(n)+ToF_GA_vector(m); 
best_coarse_time_arrival=best_coarse_time_GA+ToF_arrival_vector(t);

% EARTH
[Planet_1_kep,ksun]=uplanet(best_coarse_departure,ID1);
[ys_Planet1_coarse_departure,v_Planet1_coarse_departure]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

s_Planet1_coarse_departure = [ys_Planet1_coarse_departure';v_Planet1_coarse_departure'];
T_Planet1 = 2*pi*sqrt( Planet_1_kep(1)^3/ksun );
t_Planet1_orbit_t_span=linspace(best_coarse_departure*24*3600,best_coarse_departure*24*3600+T_Planet1,5000);

[~,ys_Planet1_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet1_orbit_t_span,s_Planet1_coarse_departure);

% SATURN
[Planet_2_kep,ksun]=uplanet(best_coarse_time_GA,ID2);
[ys_Planet2_coarse_GA,v_Planet2_coarse_GA]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

s_Planet2_coarse_GA = [ys_Planet2_coarse_GA';v_Planet2_coarse_GA'];
T_Planet2 = 2*pi*sqrt( Planet_2_kep(1)^3/ksun );
t_Planet2_orbit_t_span=linspace(best_coarse_time_GA*24*3600,best_coarse_time_GA*24*3600+T_Planet2,5000);

[~,ys_Planet2_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet2_orbit_t_span,s_Planet2_coarse_GA);

% NEO
[NEO_kep,~,~,~]=ephNEO(best_coarse_time_arrival,IDNEO);
%[NEO_kep,ksun]=uplanet(best_coarse_time_arrival,IDNEO);
[ys_NEO_coarse_arrival,v_NEO_coarse_arrival]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

s_NEO_coarse_arrival = [ys_NEO_coarse_arrival';v_NEO_coarse_arrival'];
T_NEO = 2*pi*sqrt( NEO_kep(1)^3/ksun );
t_NEO_orbit_t_span=linspace(best_coarse_time_arrival*24*3600,best_coarse_time_arrival*24*3600+T_NEO,5000);

[~,ys_NEO_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_NEO_orbit_t_span,s_NEO_coarse_arrival);

% Transfer Arc - Earth to Saturn
ToF_GA_span=linspace(best_coarse_departure*24*3600,(best_coarse_departure+ToF_GA_vector(m))*24*3600,1000);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet1_coarse_departure, ys_Planet2_coarse_GA, ToF_GA_vector(m)*24*3600, ksun, 0, 0, 0 );
s_transfer_GA=[ys_Planet1_coarse_departure',VI'];
[~,ys_trans_GA] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_GA_span,s_transfer_GA);


% Transfer Arc - Saturn to NEO
ToF_arrival_span=linspace(best_coarse_time_GA*24*3600,(best_coarse_time_GA+ToF_arrival_vector(t))*24*3600,1000);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet2_coarse_GA, ys_NEO_coarse_arrival, ToF_arrival_vector(t)*24*3600, ksun, 0, 0, 0 );
s_transfer_arrival=[ys_Planet2_coarse_GA',VI'];
[~,ys_trans_arrival] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_arrival_span,s_transfer_arrival);

% closest approach to Sun
min_dist_sun=ys_trans_arrival(:,1).^2+ys_trans_arrival(:,2).^2+ys_trans_arrival(:,3).^2; % distance squared to ease computation
min_dist_sun=sqrt(min(min_dist_sun));

fprintf('Minimum Distance to Sun: %g x Radius of Sun \n',min_dist_sun/696340);

% Plotting - HELIOCENTRIC - AU
figure(200)

% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/3;
Yt=Yt/3;
Zt=Zt/3;
image_file = 'https://upload.wikimedia.org/wikipedia/commons/9/99/Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

title('Orbit of Earth, Saturn and NEO','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color','b')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color','r')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color','g')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color','m','LineStyle','--')

% Earth to Saturn arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color','c','LineStyle','--')

legend('Earth Orbit',...
    'Saturn Orbit',...
    'NEO Orbit',...
    'Saturn Transfer',...
    'NEO transfer','Interpreter','latex','Fontsize',11)

%% Post Processing - Finding Best Solution - Fminunc

% Best solution from refinement
load("time_window_5_refinement_3.mat")
idx=find(deltaV_matrix(:)==min(deltaV_matrix,[],'all'));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

time_departure=time_departure_span(n);
ToF_GA=ToF_GA_vector(m);
ToF_arrival=ToF_arrival_vector(t);


% Best dV minimized
[best_times,deltaV_best]=fminunc(@(x) deltaV_mission(x(1),x(2),x(3)),[time_departure,ToF_GA,ToF_arrival]);

fprintf('********fminunc*****************\n')
fprintf('Best deltaV solution found in minimizer: %g\n',deltaV_best);

time_departure_minimum=mjd20002date(best_times(1));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(best_times(1)+best_times(2));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(best_times(1)+best_times(2)+best_times(3));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************** \n')
fprintf('ToF to Saturn: %g\nToF to NEO: %g \n',best_times(2),best_times(3))

[best_times,deltaV_best]=patternsearch(@(x) deltaV_mission(x(1),x(2),x(3)),[time_departure,ToF_GA,ToF_arrival]);

fprintf('********patternsearch*****************\n')
fprintf('Best deltaV solution found in minimizer: %g\n',deltaV_best);

time_departure_minimum=mjd20002date(best_times(1));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(best_times(1)+best_times(2));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(best_times(1)+best_times(2)+best_times(3));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************** \n')
fprintf('ToF to Saturn: %g\nToF to NEO: %g \n',best_times(2),best_times(3))

% [best_times,deltaV_best]=ga(@(x) deltaV_mission(x(1),x(2),x(3)),3,[-1,0,0;1,0,0;1,1,1;0,-1,0;0,1,0;0,0,-1;0,0,1],[-time_departure*0.9;time_departure*1.1;time_arrival_latest;-ToF_GA*0.9;ToF_GA*1.1;-ToF_arrival*0.9;ToF_arrival*1.1]);
% 
% fprintf('********ga*****************\n')
% fprintf('Best deltaV solution found in minimizer: %g\n',deltaV_best);
% 
% time_departure_minimum=mjd20002date(best_times(1));
% departure_disp=sprintf('%d ', time_departure_minimum);
% fprintf('Date of Departure: %s\n', departure_disp)
% 
% time_GA_minimum = mjd20002date(best_times(1)+best_times(2));
% GA_disp=sprintf('%d ', time_GA_minimum);
% fprintf('Date of GA: %s\n', GA_disp)
% 
% time_arrival_minimum = mjd20002date(best_times(1)+best_times(2)+best_times(3));
% arrival_disp=sprintf('%d ', time_arrival_minimum);
% fprintf('Date of Arrival: %s\n', arrival_disp)
% 
% fprintf('****************** \n')
% fprintf('ToF to Saturn: %g\nToF to NEO: %g \n',best_times(2),best_times(3))

[best_times,deltaV_best]=particleswarm(@(x) deltaV_mission(x(1),x(2),x(3)),3,[time_departure*0.9,ToF_GA*0.9,ToF_arrival*0.9],[time_departure*1.1,ToF_GA*1.1,ToF_arrival*1.1]);

fprintf('********particleswarm*****************\n')
fprintf('Best deltaV solution found in minimizer: %g\n',deltaV_best);

time_departure_minimum=mjd20002date(best_times(1));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(best_times(1)+best_times(2));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(best_times(1)+best_times(2)+best_times(3));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************** \n')
fprintf('ToF to Saturn: %g\nToF to NEO: %g \n',best_times(2),best_times(3))



%% Heliocentric Plot with best DeltaV result after minimization - FOR COMPARISON FIGURE - NO EXTRA PLANET PLOT


hfig = figure(200);

% CASE 1

subplot(2,2,1)

% Best solution from refinement
load("time_window_5_refinement_3.mat")
idx=find(deltaV_matrix(:)==min(deltaV_matrix,[],'all'));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

time_departure=time_departure_span(n);
ToF_GA=ToF_GA_vector(m);
ToF_arrival=ToF_arrival_vector(t);

% Plotting - HELIOCENTRIC - AU

% Best dV minimized
[best_times,deltaV_best]=fminunc(@(x) deltaV_mission(x(1),x(2),x(3)),[time_departure,ToF_GA,ToF_arrival]);
%%%% EARTH
% Plotting Orbit
[Planet_1_kep,ksun]=uplanet(best_times(1),ID1);
[ys_Planet1_best_departure,v_Planet1_best_departure]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

s_Planet1_best_departure = [ys_Planet1_best_departure';v_Planet1_best_departure'];
T_Planet1 = 2*pi*sqrt( Planet_1_kep(1)^3/ksun );
t_Planet1_orbit_t_span=linspace(best_times(1)*24*3600,best_times(1)*24*3600+T_Planet1,200);

[~,ys_Planet1_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet1_orbit_t_span,s_Planet1_best_departure);


% SATURN
[Planet_2_kep,ksun]=uplanet(best_times(1)+best_times(2),ID2);
[ys_Planet2_best_time_GA,v_Planet2_best_time_GA]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

s_Planet2_best_GA = [ys_Planet2_best_time_GA';v_Planet2_best_time_GA'];
T_Planet2 = 2*pi*sqrt( Planet_2_kep(1)^3/ksun );
t_Planet2_orbit_t_span=linspace((best_times(1)+best_times(2))*24*3600,(best_times(1)+best_times(2))*24*3600+T_Planet2,200);

[~,ys_Planet2_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet2_orbit_t_span,s_Planet2_best_GA);


% NEO
[NEO_kep,~,~,~]=ephNEO(best_times(1)+best_times(2)+best_times(3),IDNEO);
[ys_NEO_best_time_arrival,v_NEO_best_time_arrival]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

s_NEO_best_arrival = [ys_NEO_best_time_arrival';v_NEO_best_time_arrival'];
T_NEO = 2*pi*sqrt( NEO_kep(1)^3/ksun );
t_NEO_orbit_t_span=linspace((best_times(1)+best_times(2)+best_times(3))*24*3600,(best_times(1)+best_times(2)+best_times(3))*24*3600+T_NEO,200);

[~,ys_NEO_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_NEO_orbit_t_span,s_NEO_best_arrival);




% Transfer Arc - Earth to Saturn
ToF_GA_span=linspace(best_times(1)*24*3600,(best_times(1)+best_times(2))*24*3600,200);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet1_best_departure, ys_Planet2_best_time_GA, best_times(2)*24*3600, ksun, 0, 0, 0 );
s_transfer_GA=[ys_Planet1_best_departure',VI'];
[~,ys_trans_GA] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_GA_span,s_transfer_GA);


% Transfer Arc - Saturn to NEO
ToF_arrival_span=linspace((best_times(1)+best_times(2))*24*3600,(best_times(1)+best_times(2)+best_times(3))*24*3600,200);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet2_best_time_GA, ys_NEO_best_time_arrival, best_times(3)*24*3600, ksun, 0, 0, 0 );
s_transfer_arrival=[ys_Planet2_best_time_GA',VI'];
[~,ys_trans_arrival] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_arrival_span,s_transfer_arrival);

Sun_Scale = 1/3;
Saturn_Scale=1/3;
Earth_Scale=1/5;


% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Sun_Scale;
Yt=Yt*Sun_Scale;
Zt=Zt*Sun_Scale;
image_file = 'Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

%set(gca,'Color','k')

% Earth - Departure
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Earth_Scale;
Yt=Yt*Earth_Scale;
Zt=Zt*Earth_Scale;
image_file = 'Map_of_the_full_Earth.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet1_best_departure(1)/AU,Yt+ys_Planet1_best_departure(2)/AU,-Zt+ys_Planet1_best_departure(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% Saturn - Time of GA
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Saturn_Scale;
Yt=Yt*Saturn_Scale;
Zt=Zt*Saturn_Scale;
image_file = 'Map_of_the_full_Saturn.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet2_best_time_GA(1)/AU,Yt+ys_Planet2_best_time_GA(2)/AU,-Zt+ys_Planet2_best_time_GA(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% NEO  Time of Arrival
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/5;
Yt=Yt/5;
Zt=Zt/5;
image_file = 'Map_of_the_full_Mercury.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_NEO_best_time_arrival(1)/AU,Yt+ys_NEO_best_time_arrival(2)/AU,-Zt+ys_NEO_best_time_arrival(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


title('5$^{th}$ Window - Full View','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0 0.4470 0.7410],'LineStyle','--')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color',[0.4940 0.1840 0.5560],'LineStyle','-')

% Saturn to NEO arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color',[0.9290 0.6940 0.1250],'LineStyle','-')

% Limit View to Earth and NEO

% view_lim=3;
% xlim([-view_lim,view_lim])
% ylim([-view_lim,view_lim])

view([0 90])

% ZOOM CASE 1

subplot(2,2,3)

% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Sun_Scale;
Yt=Yt*Sun_Scale;
Zt=Zt*Sun_Scale;
image_file = 'Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

%set(gca,'Color','k')

% Earth - Departure
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Earth_Scale;
Yt=Yt*Earth_Scale;
Zt=Zt*Earth_Scale;
image_file = 'Map_of_the_full_Earth.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet1_best_departure(1)/AU,Yt+ys_Planet1_best_departure(2)/AU,-Zt+ys_Planet1_best_departure(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% Saturn - Time of GA
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Saturn_Scale;
Yt=Yt*Saturn_Scale;
Zt=Zt*Saturn_Scale;
image_file = 'Map_of_the_full_Saturn.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet2_best_time_GA(1)/AU,Yt+ys_Planet2_best_time_GA(2)/AU,-Zt+ys_Planet2_best_time_GA(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% NEO  Time of Arrival
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/5;
Yt=Yt/5;
Zt=Zt/5;
image_file = 'Map_of_the_full_Mercury.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_NEO_best_time_arrival(1)/AU,Yt+ys_NEO_best_time_arrival(2)/AU,-Zt+ys_NEO_best_time_arrival(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


title('5$^{th}$ Window - Zoom In','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0 0.4470 0.7410],'LineStyle','--')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color',[0.4940 0.1840 0.5560],'LineStyle','-')

% Saturn to NEO arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color',[0.9290 0.6940 0.1250],'LineStyle','-')

% Limit View to Earth and NEO

view_lim=3;
xlim([-view_lim,view_lim])
ylim([-view_lim,view_lim])

view([0 90])

% CASE 2

subplot(2,2,2)

% Best solution from refinement
load("time_window_1_refinement_3.mat")
idx=find(deltaV_matrix(:)==min(deltaV_matrix,[],'all'));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

time_departure=time_departure_span(n);
ToF_GA=ToF_GA_vector(m);
ToF_arrival=ToF_arrival_vector(t);

% Plotting - HELIOCENTRIC - AU

% Best dV minimized
[best_times,deltaV_best]=fminunc(@(x) deltaV_mission(x(1),x(2),x(3)),[time_departure,ToF_GA,ToF_arrival]);
%%%% EARTH
% Plotting Orbit
[Planet_1_kep,ksun]=uplanet(best_times(1),ID1);
[ys_Planet1_best_departure,v_Planet1_best_departure]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

s_Planet1_best_departure = [ys_Planet1_best_departure';v_Planet1_best_departure'];
T_Planet1 = 2*pi*sqrt( Planet_1_kep(1)^3/ksun );
t_Planet1_orbit_t_span=linspace(best_times(1)*24*3600,best_times(1)*24*3600+T_Planet1,200);

[~,ys_Planet1_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet1_orbit_t_span,s_Planet1_best_departure);


% SATURN
[Planet_2_kep,ksun]=uplanet(best_times(1)+best_times(2),ID2);
[ys_Planet2_best_time_GA,v_Planet2_best_time_GA]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

s_Planet2_best_GA = [ys_Planet2_best_time_GA';v_Planet2_best_time_GA'];
T_Planet2 = 2*pi*sqrt( Planet_2_kep(1)^3/ksun );
t_Planet2_orbit_t_span=linspace((best_times(1)+best_times(2))*24*3600,(best_times(1)+best_times(2))*24*3600+T_Planet2,200);

[~,ys_Planet2_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet2_orbit_t_span,s_Planet2_best_GA);


% NEO
[NEO_kep,~,~,~]=ephNEO(best_times(1)+best_times(2)+best_times(3),IDNEO);
[ys_NEO_best_time_arrival,v_NEO_best_time_arrival]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

s_NEO_best_arrival = [ys_NEO_best_time_arrival';v_NEO_best_time_arrival'];
T_NEO = 2*pi*sqrt( NEO_kep(1)^3/ksun );
t_NEO_orbit_t_span=linspace((best_times(1)+best_times(2)+best_times(3))*24*3600,(best_times(1)+best_times(2)+best_times(3))*24*3600+T_NEO,200);

[~,ys_NEO_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_NEO_orbit_t_span,s_NEO_best_arrival);


% Transfer Arc - Earth to Saturn
ToF_GA_span=linspace(best_times(1)*24*3600,(best_times(1)+best_times(2))*24*3600,100);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet1_best_departure, ys_Planet2_best_time_GA, best_times(2)*24*3600, ksun, 0, 0, 0 );
s_transfer_GA=[ys_Planet1_best_departure',VI'];
[~,ys_trans_GA] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_GA_span,s_transfer_GA);


% Transfer Arc - Saturn to NEO
ToF_arrival_span=linspace((best_times(1)+best_times(2))*24*3600,(best_times(1)+best_times(2)+best_times(3))*24*3600,100);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet2_best_time_GA, ys_NEO_best_time_arrival, best_times(3)*24*3600, ksun, 0, 0, 0 );
s_transfer_arrival=[ys_Planet2_best_time_GA',VI'];
[~,ys_trans_arrival] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_arrival_span,s_transfer_arrival);

Sun_Scale = 1/3;
Saturn_Scale=1/3;
Earth_Scale=1/5;


% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Sun_Scale;
Yt=Yt*Sun_Scale;
Zt=Zt*Sun_Scale;
image_file = 'Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

%set(gca,'Color','k')

% Earth - Departure
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Earth_Scale;
Yt=Yt*Earth_Scale;
Zt=Zt*Earth_Scale;
image_file = 'Map_of_the_full_Earth.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet1_best_departure(1)/AU,Yt+ys_Planet1_best_departure(2)/AU,-Zt+ys_Planet1_best_departure(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% Saturn - Time of GA
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Saturn_Scale;
Yt=Yt*Saturn_Scale;
Zt=Zt*Saturn_Scale;
image_file = 'Map_of_the_full_Saturn.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet2_best_time_GA(1)/AU,Yt+ys_Planet2_best_time_GA(2)/AU,-Zt+ys_Planet2_best_time_GA(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% NEO  Time of Arrival
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/5;
Yt=Yt/5;
Zt=Zt/5;
image_file = 'Map_of_the_full_Mercury.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_NEO_best_time_arrival(1)/AU,Yt+ys_NEO_best_time_arrival(2)/AU,-Zt+ys_NEO_best_time_arrival(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


title('1$^{st}$ Window - Full View','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0 0.4470 0.7410],'LineStyle','--')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color',[0.4940 0.1840 0.5560],'LineStyle','-')

% Saturn to NEO arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color',[0.9290 0.6940 0.1250],'LineStyle','-')

% Limit View to Earth and NEO

% view_lim=3;
% xlim([-view_lim,view_lim])
% ylim([-view_lim,view_lim])

view([0 90])

% ZOOM Case 2

subplot(2,2,4)

% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Sun_Scale;
Yt=Yt*Sun_Scale;
Zt=Zt*Sun_Scale;
image_file = 'Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

%set(gca,'Color','k')

% Earth - Departure
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Earth_Scale;
Yt=Yt*Earth_Scale;
Zt=Zt*Earth_Scale;
image_file = 'Map_of_the_full_Earth.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet1_best_departure(1)/AU,Yt+ys_Planet1_best_departure(2)/AU,-Zt+ys_Planet1_best_departure(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% Saturn - Time of GA
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Saturn_Scale;
Yt=Yt*Saturn_Scale;
Zt=Zt*Saturn_Scale;
image_file = 'Map_of_the_full_Saturn.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet2_best_time_GA(1)/AU,Yt+ys_Planet2_best_time_GA(2)/AU,-Zt+ys_Planet2_best_time_GA(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


% NEO  Time of Arrival
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/5;
Yt=Yt/5;
Zt=Zt/5;
image_file = 'Map_of_the_full_Mercury.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_NEO_best_time_arrival(1)/AU,Yt+ys_NEO_best_time_arrival(2)/AU,-Zt+ys_NEO_best_time_arrival(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


title('1$^{st}$ Window - Zoom In','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0 0.4470 0.7410],'LineStyle','--')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color',[0.4940 0.1840 0.5560],'LineStyle','-')

% Saturn to NEO arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color',[0.9290 0.6940 0.1250],'LineStyle','-')

% Limit View to Earth and NEO

view_lim=3;
xlim([-view_lim,view_lim])
ylim([-view_lim,view_lim])

view([0 90])

sgtitle('Comparison between minimized 1$^{st}$ and 5$^{th}$ departure window','Interpreter','latex','Fontsize',11)

% legend('Earth Orbit',...
%     'Saturn Orbit',...
%     'NEO Orbit',...
%     'Saturn Transfer',...
%     'NEO transfer','Interpreter','latex','Fontsize',11)

fname = 'Helio_minimized_comparison';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',15) % adjust fontsize to your document

%set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')

%% Heliocentric Plot with best DeltaV result after minimization

%%%% EARTH
% Plotting Orbit
[Planet_1_kep,ksun]=uplanet(best_times(1),ID1);
[ys_Planet1_best_departure,v_Planet1_best_departure]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

s_Planet1_best_departure = [ys_Planet1_best_departure';v_Planet1_best_departure'];
T_Planet1 = 2*pi*sqrt( Planet_1_kep(1)^3/ksun );
t_Planet1_orbit_t_span=linspace(best_times(1)*24*3600,best_times(1)*24*3600+T_Planet1,200);

[~,ys_Planet1_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet1_orbit_t_span,s_Planet1_best_departure);

% Plotting position with surface    
% Earth at Time of Gravity Assist
[Planet_1_kep,ksun]=uplanet(best_times(1)+best_times(2),ID1);
[ys_Planet1_best_time_GA,~]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);

% Earth at Time of NEO arrival
[Planet_1_kep,ksun]=uplanet(best_times(1)+best_times(2)+best_times(3),ID1);
[ys_Planet1_best_time_arrival,~]=kep2car(Planet_1_kep(1),...
    Planet_1_kep(2),...
    rad2deg(Planet_1_kep(3)),...
    rad2deg(Planet_1_kep(4)),...
    rad2deg(Planet_1_kep(5)),...
    rad2deg(Planet_1_kep(6)),...
    ksun);


% SATURN
[Planet_2_kep,ksun]=uplanet(best_times(1)+best_times(2),ID2);
[ys_Planet2_best_time_GA,v_Planet2_best_time_GA]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

s_Planet2_best_GA = [ys_Planet2_best_time_GA';v_Planet2_best_time_GA'];
T_Planet2 = 2*pi*sqrt( Planet_2_kep(1)^3/ksun );
t_Planet2_orbit_t_span=linspace((best_times(1)+best_times(2))*24*3600,(best_times(1)+best_times(2))*24*3600+T_Planet2,200);

[~,ys_Planet2_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_Planet2_orbit_t_span,s_Planet2_best_GA);

% Plotting position with surface    
% Saturn at Time of Departure
[Planet_2_kep,~]=uplanet(best_times(1),ID2);
[ys_Planet2_best_departure,~]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

% Saturn at Time of Arrival
[Planet_2_kep,~]=uplanet(best_times(1)+best_times(2)+best_times(3),ID2);
[ys_Planet2_best_time_arrival,~]=kep2car(Planet_2_kep(1),...
    Planet_2_kep(2),...
    rad2deg(Planet_2_kep(3)),...
    rad2deg(Planet_2_kep(4)),...
    rad2deg(Planet_2_kep(5)),...
    rad2deg(Planet_2_kep(6)),...
    ksun);

% NEO
[NEO_kep,~,~,~]=ephNEO(best_times(1)+best_times(2)+best_times(3),IDNEO);
[ys_NEO_best_time_arrival,v_NEO_best_time_arrival]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

s_NEO_best_arrival = [ys_NEO_best_time_arrival';v_NEO_best_time_arrival'];
T_NEO = 2*pi*sqrt( NEO_kep(1)^3/ksun );
t_NEO_orbit_t_span=linspace((best_times(1)+best_times(2)+best_times(3))*24*3600,(best_times(1)+best_times(2)+best_times(3))*24*3600+T_NEO,200);

[~,ys_NEO_orbit_complete] = ode89(@(t,y) ode_2bp(t,y,ksun),t_NEO_orbit_t_span,s_NEO_best_arrival);

% Plotting position with surface    
% NEO at Time of Departure
[NEO_kep,~,~,~]=ephNEO(best_times(1),IDNEO);
[ys_NEO_best_departure,~]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);

% NEO at Time of GA
[NEO_kep,~,~,~]=ephNEO(best_times(1)+best_times(2),IDNEO);
[ys_NEO_best_time_GA,~]=kep2car(NEO_kep(1),...
    NEO_kep(2),...
    rad2deg(NEO_kep(3)),...
    rad2deg(NEO_kep(4)),...
    rad2deg(NEO_kep(5)),...
    rad2deg(NEO_kep(6)),...
    ksun);


% Transfer Arc - Earth to Saturn
ToF_GA_span=linspace(best_times(1)*24*3600,(best_times(1)+best_times(2))*24*3600,100);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet1_best_departure, ys_Planet2_best_time_GA, best_times(2)*24*3600, ksun, 0, 0, 0 );
s_transfer_GA=[ys_Planet1_best_departure',VI'];
[ts_trans_GA,ys_trans_GA] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_GA_span,s_transfer_GA);


% Transfer Arc - Saturn to NEO
ToF_arrival_span=linspace((best_times(1)+best_times(2))*24*3600,(best_times(1)+best_times(2)+best_times(3))*24*3600,100);
[~,~,~,~,VI,~,~,~] = lambertMR( ys_Planet2_best_time_GA, ys_NEO_best_time_arrival, best_times(3)*24*3600, ksun, 0, 0, 0 );
s_transfer_arrival=[ys_Planet2_best_time_GA',VI'];
[ts_trans_arrival,ys_trans_arrival] = ode89(@(t,y) ode_2bp(t,y,ksun),ToF_arrival_span,s_transfer_arrival);

% Plotting - HELIOCENTRIC - AU
hfig=figure(200);

Sun_Scale = 1/3;
Saturn_Scale=1/3;
Earth_Scale=1/5;


% SUN
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Sun_Scale;
Yt=Yt*Sun_Scale;
Zt=Zt*Sun_Scale;
image_file = 'Map_of_the_full_sun.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

%set(gca,'Color','k')

% Earth - 3 positions
% Departure
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Earth_Scale;
Yt=Yt*Earth_Scale;
Zt=Zt*Earth_Scale;
image_file = 'Map_of_the_full_Earth.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet1_best_departure(1)/AU,Yt+ys_Planet1_best_departure(2)/AU,-Zt+ys_Planet1_best_departure(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Time of GA
Globe = surf(Xt+ys_Planet1_best_time_GA(1)/AU,Yt+ys_Planet1_best_time_GA(2)/AU,-Zt+ys_Planet1_best_time_GA(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Time of Arrival
Globe = surf(Xt+ys_Planet1_best_time_arrival(1)/AU,Yt+ys_Planet1_best_time_arrival(2)/AU,-Zt+ys_Planet1_best_time_arrival(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Saturn - 3 positions
% Departure
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt*Saturn_Scale;
Yt=Yt*Saturn_Scale;
Zt=Zt*Saturn_Scale;
image_file = 'Map_of_the_full_Saturn.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_Planet2_best_departure(1)/AU,Yt+ys_Planet2_best_departure(2)/AU,-Zt+ys_Planet2_best_departure(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Time of GA
Globe = surf(Xt+ys_Planet2_best_time_GA(1)/AU,Yt+ys_Planet2_best_time_GA(2)/AU,-Zt+ys_Planet2_best_time_GA(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Time of Arrival
Globe = surf(Xt+ys_Planet2_best_time_arrival(1)/AU,Yt+ys_Planet2_best_time_arrival(2)/AU,-Zt+ys_Planet2_best_time_arrival(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% NEO - 3 positions
% Departure
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks
Xt=Xt/5;
Yt=Yt/5;
Zt=Zt/5;
image_file = 'Map_of_the_full_Mercury.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt+ys_NEO_best_departure(1)/AU,Yt+ys_NEO_best_departure(2)/AU,-Zt+ys_NEO_best_departure(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Time of GA
Globe = surf(Xt+ys_NEO_best_time_GA(1)/AU,Yt+ys_NEO_best_time_GA(2)/AU,-Zt+ys_NEO_best_time_GA(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Time of Arrival
Globe = surf(Xt+ys_NEO_best_time_arrival(1)/AU,Yt+ys_NEO_best_time_arrival(2)/AU,-Zt+ys_NEO_best_time_arrival(3)/AU); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';


title('Best Mission to NEO','Interpreter','latex','Fontsize',11);
xlabel('$x$ [AU]','Interpreter','latex','Fontsize',11); 
ylabel('$y$ [AU]','Interpreter','latex','Fontsize',11); 
zlabel('$z$ [AU]','Interpreter','latex','Fontsize',11);
axis equal;
grid on;

% Earth Orbit
plot3(ys_Planet1_orbit_complete(:,1)/AU,ys_Planet1_orbit_complete(:,2)/AU,ys_Planet1_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0 0.4470 0.7410],'LineStyle','--')

% Saturn Orbit
plot3(ys_Planet2_orbit_complete(:,1)/AU,ys_Planet2_orbit_complete(:,2)/AU,ys_Planet2_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')

% NEO Orbit
plot3(ys_NEO_orbit_complete(:,1)/AU,ys_NEO_orbit_complete(:,2)/AU,ys_NEO_orbit_complete(:,3)/AU,'LineWidth',2,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')

% Earth to Saturn arc
plot3(ys_trans_GA(:,1)/AU,ys_trans_GA(:,2)/AU,ys_trans_GA(:,3)/AU,'LineWidth',2,'Color',[0.4940 0.1840 0.5560],'LineStyle','-')

% Saturn to NEO arc
plot3(ys_trans_arrival(:,1)/AU,ys_trans_arrival(:,2)/AU,ys_trans_arrival(:,3)/AU,'LineWidth',2,'Color',[0.9290 0.6940 0.1250],'LineStyle','-')

% Annotations
% Earth
text(ys_Planet1_best_departure(1)/AU+0.5,ys_Planet1_best_departure(2)/AU,ys_Planet1_best_departure(3)/AU,'Dep','Color','blue','FontSize',10)
text(ys_Planet1_best_time_GA(1)/AU+0.5,ys_Planet1_best_time_GA(2)/AU,ys_Planet1_best_time_GA(3)/AU,'GA','Color','blue','FontSize',10)
text(ys_Planet1_best_time_arrival(1)/AU+0.5,ys_Planet1_best_time_arrival(2)/AU,ys_Planet1_best_time_arrival(3)/AU,'Arr','Color','blue','FontSize',10)

% Saturn
text(ys_Planet2_best_departure(1)/AU+0.8,ys_Planet2_best_departure(2)/AU,ys_Planet2_best_departure(3)/AU,'   Dep','Color','red','FontSize',10)
text(ys_Planet2_best_time_GA(1)/AU-0.2,ys_Planet2_best_time_GA(2)/AU-1,ys_Planet2_best_time_GA(3)/AU,'   GA','Color','red','FontSize',10)
text(ys_Planet2_best_time_arrival(1)/AU+0.7,ys_Planet2_best_time_arrival(2)/AU,ys_Planet2_best_time_arrival(3)/AU,'   Arr','Color','red','FontSize',10)

% NEO
text(ys_NEO_best_departure(1)/AU-0.4,ys_NEO_best_departure(2)/AU+0.6,ys_NEO_best_departure(3)/AU,'Dep','Color',[1 50 32]/258,'FontSize',10)
text(ys_NEO_best_time_GA(1)/AU,ys_NEO_best_time_GA(2)/AU,ys_NEO_best_time_GA(3)/AU,' GA','Color',[1 50 32]/258,'FontSize',10)
text(ys_NEO_best_time_arrival(1)/AU-0.4,ys_NEO_best_time_arrival(2)/AU-0.6,ys_NEO_best_time_arrival(3)/AU,'Arr','Color',[1 50 32]/258,'FontSize',10)

% Limit View to Earth and NEO

% view_lim=3;
% xlim([-view_lim,view_lim])
% ylim([-view_lim,view_lim])

ylim([-2,10])

view([0 90])

legend('Earth Orbit',...
    'Saturn Orbit',...
    'NEO Orbit',...
    'Saturn Transfer',...
    'NEO transfer','Interpreter','latex','Fontsize',11)

fname = 'Best_Helio';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',15) % adjust fontsize to your document

%set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')

%% Mission Parameters

R_saturn=astroConstants(26);
mu_Saturn=astroConstants(16);
rp_crit=1.3*R_saturn;

time_GA = best_times(1)+best_times(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position and Velocity of Saturn
[Saturn_kep,~]=uplanet(time_GA,ID2);

[~,V_Saturn_tGA]=kep2car(Saturn_kep(1),...
    Saturn_kep(2),...
    rad2deg(Saturn_kep(3)),...
    rad2deg(Saturn_kep(4)),...
    rad2deg(Saturn_kep(5)),...
    rad2deg(Saturn_kep(6)),...
    ksun);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,deltaV_departure,~,~,V_inf_minus]=delta_V_Interplanetary(time_departure,time_GA,ID1,ID2);

time_arrival    = best_times(1)+best_times(2)+best_times(3);

[~,~,V_inf_plus,deltaV_arrival,~]=delta_V_Interplanetary(time_GA,time_arrival,ID2,IDNEO);

v_inf_minus = V_inf_minus-V_Saturn_tGA;
v_inf_plus  = V_inf_plus-V_Saturn_tGA;


delta_angle = acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_minus)*norm(v_inf_plus)));
rp_specific = @(rp) rp_solver(rp,norm(v_inf_minus),norm(v_inf_plus),delta_angle,mu_Saturn);

[rp,~,rp_sol_flag,~] = fzero(rp_specific,R_saturn*10,optimset('Display','off'));

v_inf_minus_sq_norm = sum(v_inf_minus.*v_inf_minus);
v_inf_plus_sq_norm  = sum(v_inf_plus.*v_inf_plus);

vp_minus_norm     = sqrt(((rp*v_inf_minus_sq_norm+2*mu_Saturn)/rp));
vp_plus_norm      = sqrt(((rp*v_inf_plus_sq_norm+2*mu_Saturn)/rp));

ecc_minus=(mu_Saturn+rp*norm(v_inf_minus)^2)/mu_Saturn;
ecc_plus=(mu_Saturn+rp*norm(v_inf_plus)^2)/mu_Saturn;

delta_minus=2*asin(1/ecc_minus);
delta_plus=2*asin(1/ecc_plus);

a_minus=rp/(1-ecc_minus);
a_plus=rp/(1-ecc_plus);

impact_parameter=-a_minus*ecc_minus*cos(delta_angle/2);

deltaV_flyby=norm(v_inf_plus-v_inf_minus);

deltaV_GA   = abs(vp_plus_norm-vp_minus_norm);

deltaV_total=deltaV_GA+norm(deltaV_departure)+norm(deltaV_arrival);

fprintf('Best Mission Parameters:\n');

time_departure_minimum=mjd20002date(best_times(1));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(best_times(1)+best_times(2));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(best_times(1)+best_times(2)+best_times(3));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************** \n')
fprintf('ToF to Saturn: %g\nToF to NEO: %g \n',best_times(2),best_times(3))

fprintf('****************** \n')
fprintf('Gravity Assist:\n')
fprintf('Eccentricity minus: %g\nEccentricity plus: %g\n',ecc_minus,ecc_plus)
fprintf('Semi-major Axis minus: %g\nSemi-major Axis plus: %g\n',a_minus,a_plus)
fprintf('Delta Axis minus: %g\nDelta plus: %g\n',rad2deg(delta_minus),rad2deg(delta_plus))
fprintf('Impact Parameter:%g x Radius of Saturn\n',impact_parameter/R_saturn)
fprintf('Radius of Perigee: %g [km]\nRadius of Perigee: %g x Radius of Saturn\n',rp,rp/R_saturn)
fprintf('DeltaV of GA: %g [km/s]\n',deltaV_GA)
fprintf('DeltaV of flyby: %g [km/s]\n',deltaV_flyby)
fprintf('DeltaV of mission: %g [km/s]\n',deltaV_total)
fprintf('DeltaV of departure: %g [km/s]\n',norm(deltaV_departure))
fprintf('DeltaV of arrival: %g [km/s]\n',norm(deltaV_arrival))

%% Saturn Fly by Plot

% Vectors expressed in the heliocentric ecliptic frame

% Curve Direction. Plane of Hyperbola in Heliocentric
% Given by the normalized cross product of v_inf_minus and v_inf_plus
u   = cross(v_inf_minus,v_inf_plus);
u   = u/norm(u);

theta_inf_minus   = pi-acos(1/ecc_minus); % Value must be greater than 90º    
beta_angle_minus  = pi-theta_inf_minus; 

theta_inf_plus    = pi-acos(1/ecc_plus);  % Value must be greater than 90º    
beta_angle_plus   = pi-theta_inf_plus; 

% Hyperbola in Perifocal Frame

% Planetocentric Hyperbola is defined in Perifocal frame:
%   x - Pericenter
%   z - alligned with Angular momentum // Plane 
%   y - Orthogonal

% Must solve ODE backwards and forward in time starting from pericenter
figure(500)
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks 
image_file = 'Map_of_the_full_Saturn.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

title('Perifocal Saturn Centered Powered Fly-By','Interpreter','latex','Fontsize',11);
xlabel('$x [R_{\oplus}]$','Interpreter','latex','Fontsize',11); 
ylabel('$y [R_{\oplus}]$','Interpreter','latex','Fontsize',11); 
zlabel('$z [R_{\oplus}]$','Interpreter','latex','Fontsize',11);
axis equal;
grid on;
hold on

vp_minus    = vp_minus_norm*[0;1;0];   % Velocity at Perigee is aligned from y+ axis
vp_plus     = vp_plus_norm*[0;1;0];

s_incoming  = [rp*[1;0;0];vp_minus];
s_outgoing  = [rp*[1;0;0];vp_plus];

% Solve backwards and forward in time
precision=1000;
t_span=linspace(0,200000,precision);

[~,ys_back] = ode89(@(t,y) ode_2bp(t,y,mu_Saturn),-t_span,s_incoming);
%ys_back=flipud(ys_back); % flip so last point is perigee, not backwards
[~,ys_front] = ode89(@(t,y) ode_2bp(t,y,mu_Saturn),t_span,s_outgoing);

plot3(ys_back(:,1)/R_saturn,ys_back(:,2)/R_saturn,ys_back(:,3)/R_saturn,'LineWidth',2)
plot3(ys_front(:,1)/R_saturn,ys_front(:,2)/R_saturn,ys_front(:,3)/R_saturn,'LineWidth',2)

% Asymptotes
% Confluence point if two asymptotes C in Perifocal is rp+a. Create a line
% defined by this point and a vector parallel to v_inf_minus and second
% asymptote is same procedure, with other C and v_inf_plus

C_minus = (rp+abs(a_minus))/R_saturn;  % Confluence point between apse line and asymptote - Incoming
C_plus  = (rp+abs(a_plus))/R_saturn;   % Confluence point between apse line and asymptote - Outgoing

asymptote_minus = [C_minus*[1;0;0],C_minus*[1;0;0]+[cos(-theta_inf_minus);sin(-theta_inf_minus);0]*(40)];
asymptote_plus = [C_plus*[1;0;0],C_plus*[1;0;0]+[cos(theta_inf_plus);sin(theta_inf_plus);0]*(40)];

% Plot asymptotes
plot3(asymptote_minus(1,:),asymptote_minus(2,:),asymptote_minus(3,:),'LineWidth',1.5)
plot3(asymptote_plus(1,:),asymptote_plus(2,:),asymptote_plus(3,:),'LineWidth',1.5)

plot3([-10,C_minus],[0,0],[0,0],'--','LineWidth',1.5,'Color','k')
% xlim([-6,6])
% ylim([-9,9])

legend('Incoming Hyperbola',...
    'Outgoing Hyperbola',...
    'Asymptote for Incoming hyperbola',...
    'Asymptote for Outgoing hyperbola',...
    '','Interpreter','latex','Fontsize',11)

view(0,90)

%% Hyperbola in Heliocentric Frame

% Plot Planet
hfig=figure(501);
[Xt,Yt,Zt] = sphere(100); %Creates an sphere wit 100 spherical blocks 
image_file = 'Map_of_the_full_Saturn.jpg';
[atlante] = imread(image_file); %Used to print the image in the surface 'earth_globe.JPG'
Globe = surf(Xt,Yt,-Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
Globe.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

title('Powered Gravity Assist','Interpreter','latex','Fontsize',11);
xlabel('$x [R_{\oplus}]$','Interpreter','latex','Fontsize',11); 
ylabel('$y [R_{\oplus}]$','Interpreter','latex','Fontsize',11); 
zlabel('$z [R_{\oplus}]$','Interpreter','latex','Fontsize',11);
axis equal;
grid on;
hold on

% Now must rotate rp and vp to be in heliocentric frame. 
% Perifocal X axis is coincident with Rp, which is inclined Beta radians
% wrt to v_inf_minus. So, create a vector in the direction of
% v_inf_minus and magnitud of rp. Then pass from v_inf_minus direction to
% heliocentric rp by rotating this vector -Beta around plane u.

rp_helio    = v_inf_minus/vecnorm(v_inf_minus)*vecnorm(rp); % rp in direction of v_inf_minus
rp_helio    = v_rotate_rodriguez(rp_helio,u,-beta_angle_minus); % perigee rotated to heliocentric frame

% Velocity at perigee is always at 90º wrt rp_helio, equivalent to
% direction of v_inf_minus rotated -Beta+Pi/2

% Incoming Parabola
vp_minus_helio    = v_inf_minus/vecnorm(v_inf_minus)*vp_minus_norm;
vp_minus_helio    = v_rotate_rodriguez(vp_minus_helio,u,-beta_angle_minus+pi/2);

s_minus=[rp_helio';vp_minus_helio'];

% Outgoing Hyperbola
vp_plus_helio     = v_inf_minus/vecnorm(v_inf_minus)*vp_plus_norm;
vp_plus_helio     = v_rotate_rodriguez(vp_plus_helio,u,-beta_angle_minus+pi/2);
% Procedure can be done the same, just change magnitude

s_plus=[rp_helio';vp_plus_helio'];

% Solve backwards and forward in time
precision=20;
t_span=linspace(0,200000,precision);


[~,ys_back] = ode89(@(t,y) ode_2bp(t,y,mu_Saturn),-t_span,s_minus);
%ys_back=flipud(ys_back); % flip so last point is perigee, not backwards
[~,ys_front] = ode89(@(t,y) ode_2bp(t,y,mu_Saturn),t_span,s_plus);

% Hyperbola
plot3(ys_back(:,1)/R_saturn,ys_back(:,2)/R_saturn,ys_back(:,3)/R_saturn,'LineWidth',2)
plot3(ys_front(:,1)/R_saturn,ys_front(:,2)/R_saturn,ys_front(:,3)/R_saturn,'LineWidth',2)

% Asymptotes
C_minus     = (rp_helio/norm(rp_helio))*((rp+abs(a_minus))/R_saturn); % Unit Vector * Length
C_plus      = (rp_helio/norm(rp_helio))*((rp+abs(a_plus))/R_saturn);

asymptote_minus = [C_minus;C_minus+(v_inf_minus/norm(v_inf_minus))*(-40)];
%asymptote_minus = asymptote_minus';

asymptote_plus  = [C_plus;C_plus+(v_inf_plus/norm(v_inf_plus))*(+40)];
%asymptote_plus  = asymptote_plus';

% Plot Asymptotes
plot3(asymptote_minus(:,1),asymptote_minus(:,2),asymptote_minus(:,3),'LineWidth',1.5)
plot3(asymptote_plus(:,1),asymptote_plus(:,2),asymptote_plus(:,3),'LineWidth',1.5)

% Apse line
apse_2=C_minus+(rp_helio/norm(rp_helio))*(-30);

plot3([C_minus(1),apse_2(1)],[C_minus(2),apse_2(2)],[C_minus(3),apse_2(3)],'--','LineWidth',1.5,'Color','k')

% scatter3(rp_helio(1)/R_saturn,rp_helio(2)/R_saturn,rp_helio(3)/R_saturn,40)

% Obit Plane - Debugging Purposes. All lines should be contained in this
% plane

%   w = null(u); % Find two orthonormal vectors which are orthogonal to v
%    [P,Q] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
%    X = (C_minus(1))+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
%    Y = (C_minus(2))+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
%    Z = (C_minus(3))+w(3,1)*P+w(3,2)*Q;
%    surf(X,Y,Z)

% Saturn Velocity
[Saturn_kep,~]=uplanet(time_GA,ID2);

[~,V_Saturn_tGA]=kep2car(Saturn_kep(1),...
    Saturn_kep(2),...
    rad2deg(Saturn_kep(3)),...
    rad2deg(Saturn_kep(4)),...
    rad2deg(Saturn_kep(5)),...
    rad2deg(Saturn_kep(6)),...
    ksun);

mArrow3([0 0 0],V_Saturn_tGA); % Saturn Velocity Vector

% legend('Incoming Hyperbola',...
%     'Outgoing Hyperbola',...
%     'Asymptote for Incoming hyperbola',...
%     'Asymptote for Outgoing hyperbola',...
%     '','','Interpreter','latex','Fontsize',11)

legend('Outgoing Hyperbola',...
    'Incoming Hyperbola',...
    'Asymptote for Outgoing hyperbola',...
    'Asymptote for Incoming hyperbola',...
    '','','Interpreter','latex','Fontsize',11)


view(0,90)

fname = 'Hyperbola_Helio';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',15) % adjust fontsize to your document

%set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')


%% Post Processing - 2nd Refinement - Best time windows 4D plot

load('time_window_3_long_refinement_2.mat')
best_limit=1000;
cmap=colormap(jet(best_limit));

best_deltaV=sort(deltaV_matrix(:)); % Best deltaV
best_deltaV=best_deltaV(1:best_limit);

[~,idx_array]=sort(deltaV_matrix(:)); % CREATES A SORTED COLUMN VECTOR OF EXISTING DELTAV values

[n2, m2, t2] = ind2sub(size(deltaV_matrix),idx_array);

n2=n2(1:best_limit);
m2=m2(1:best_limit);
t2=t2(1:best_limit);


[minDeltaV, idx] = min(deltaV_matrix(:));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

% min DeltaV
fprintf('Minimum deltaV: %g \n',minDeltaV)

time_departure_minimum=mjd20002date(time_departure_span(n));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m)+ToF_arrival_vector(t));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************************** \n');
fprintf('Radius of periguee: %g x Radius of Saturn\nTurn Angle (deg): %g\nToF to Saturn (days): %g\nToF to NEO (days): %g\n',rp_matrix(n,m,t)/R_saturn,rad2deg(delta_angle_matrix(n,m,t)),ToF_GA_vector(m),ToF_arrival_vector(t))


hfig=figure(402);
scatter3(time_departure_span(n2),ToF_GA_vector(m2),ToF_arrival_vector(t2),2,cmap,'filled','MarkerFaceAlpha',.2)
clim([16,20])
colorbar
xlim([9900,12100])
hold on



for k=1:6
    load_string=['time_window_',num2str(k),'_refinement_2.mat'];
    load(load_string);


    best_deltaV=sort(deltaV_matrix(:)); % Best deltaV
    best_deltaV=best_deltaV(1:best_limit);

    [~,idx_array]=sort(deltaV_matrix(:)); % CREATES A SORTED COLUMN VECTOR OF EXISTING DELTAV values

    [n2, m2, t2] = ind2sub(size(deltaV_matrix),idx_array);

    n2=n2(1:best_limit);
    m2=m2(1:best_limit);
    t2=t2(1:best_limit);

    scatter3(time_departure_span(n2),ToF_GA_vector(m2),ToF_arrival_vector(t2),2,cmap,'filled','MarkerFaceAlpha',.2)
    clim([16,20])
    colorbar


[minDeltaV, idx] = min(deltaV_matrix(:));
[n, m, t] = ind2sub(size(deltaV_matrix),idx);

% min DeltaV
fprintf('---------------------------\n');
fprintf('Time Window number: %g\n',k)

fprintf('Minimum deltaV: %g \n',minDeltaV);

time_departure_minimum=mjd20002date(time_departure_span(n));
departure_disp=sprintf('%d ', time_departure_minimum);
fprintf('Date of Departure: %s\n', departure_disp)

time_GA_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m));
GA_disp=sprintf('%d ', time_GA_minimum);
fprintf('Date of GA: %s\n', GA_disp)

time_arrival_minimum = mjd20002date(time_departure_span(n)+ToF_GA_vector(m)+ToF_arrival_vector(t));
arrival_disp=sprintf('%d ', time_arrival_minimum);
fprintf('Date of Arrival: %s\n', arrival_disp)

fprintf('****************************** \n');
fprintf('ToF to Saturn (days): %g\nToF to NEO (days): %g\n',ToF_GA_vector(m),ToF_arrival_vector(t))


end

title('Mission $\Delta v$ - Refined - Complete time span','Interpreter','latex','Fontsize',11);
xlabel('departure time [days - mjd2000]','Interpreter','latex','Fontsize',11); 
ylabel('Time of Flight to Saturn','Interpreter','latex','Fontsize',11); 
zlabel('Time of Flight to NEO','Interpreter','latex','Fontsize',11);
%xlim([time_departure_span(1),time_departure_span(end)])
%xlim([time_departure_earliest,time_departure_latest])
% axis equal;
grid on;

% fname = 'Refined_Grid';
% 
% picturewidth = 20; % set this parameter and keep it forever
% hw_ratio = 0.65; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% 
% %set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hfig,fname,'-dpdf','-painters','-fillpage')
% %print(hfig,fname,'-dpng','-painters')

%% Minimizing Function

function deltaV_total=deltaV_mission(time_departure,ToF_GA,ToF_arrival)

ID1=3;
ID2=6;
IDNEO=51;

ksun=astroConstants(4);

R_saturn=astroConstants(26);
mu_Saturn=astroConstants(16);
rp_crit=1.3*R_saturn;

time_GA = time_departure+ToF_GA;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Position and Velocity of Saturn
        [Saturn_kep,~]=uplanet(time_GA,ID2);

        [~,V_Saturn_tGA]=kep2car(Saturn_kep(1),...
            Saturn_kep(2),...
            rad2deg(Saturn_kep(3)),...
            rad2deg(Saturn_kep(4)),...
            rad2deg(Saturn_kep(5)),...
            rad2deg(Saturn_kep(6)),...
            ksun);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,deltaV_departure,~,~,V_inf_minus]=delta_V_Interplanetary(time_departure,time_GA,ID1,ID2);

time_arrival    = time_GA+ToF_arrival; 

[~,~,V_inf_plus,deltaV_arrival,~]=delta_V_Interplanetary(time_GA,time_arrival,ID2,IDNEO);

v_inf_minus = V_inf_minus-V_Saturn_tGA;     
v_inf_plus  = V_inf_plus-V_Saturn_tGA;

delta_angle = acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_minus)*norm(v_inf_plus)));     
rp_specific = @(rp) rp_solver(rp,norm(v_inf_minus),norm(v_inf_plus),delta_angle,mu_Saturn);

[rp,~,rp_sol_flag,~] = fzero(rp_specific,R_saturn*10,optimset('Display','off'));

if (rp<rp_crit) || (rp_sol_flag<1)
    warning('Rp solver not correct - Error flag or rp<rpcrit')
end

v_inf_minus_sq_norm = sum(v_inf_minus.*v_inf_minus); 
v_inf_plus_sq_norm  = sum(v_inf_plus.*v_inf_plus);

vp_minus_norm     = sqrt(((rp*v_inf_minus_sq_norm+2*mu_Saturn)/rp)); 
vp_plus_norm      = sqrt(((rp*v_inf_plus_sq_norm+2*mu_Saturn)/rp));  

deltaV_GA   = abs(vp_plus_norm-vp_minus_norm);

deltaV_total=deltaV_GA+norm(deltaV_departure)+norm(deltaV_arrival);

end