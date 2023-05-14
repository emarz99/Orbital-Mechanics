function[v_inf_i, v_inf_f, DV_cost_vec, Time_in_SOI, delta_deg, hp]=Hyperbola_plot(V_minus, V_plus, fb_Planet, gaTime)

% Hyperbola_plot.m - plot in 3Dspace of the incoming and outcoming
%                    hyperbolas attached at their pericenyter
%
% PROTOTYPE:
%    [v_inf_i, v_inf_f, DV_cost_vec, Time_in_SOI, hp, delta_deg]=Hyperbola_plot(V_minus, V_plus, fb_Planet, gaTime)
%
% DESCRIPTION: taking as input gravity assist Planet, date of maneuver and
%              Velocity of S/C on the interplanetary leg before and after fly-by define
%              all the fly-by parameters and plots it
%
% INPUT:        
%   V_minus[3]     incoming velocity of S/C in heliocentric ref frame 
%   V_plus[3]      outcoming velocity of S/C in heliocentric ref frame 
%   fb_planet[1]             Integer number identifying the GA body 
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
%   gaTime[6]      Date in the Gregorian calendar, as a 6-element vector
%                       [year, month, day, hour, minute, second]
%   
%
%
% OUTPUT:
%  v_inf_i[3]      incoming velocity vector of S/C wrt to planet obtain as 
%                               V_minus-Vplanet
%  v_inf_i[3]      outcoming velocity vector of S/C wrt to planet obtain as 
%                               V_plus-Vplanet
%  DV_cost_vec[3]  cost of gravity assist obtained represented to the DV at
%                  the pericenter to match the  pericenter velocity of the 
%                  incoming and outcoming arc
%  Time_in_SOI[1]  time spent inside SOI of selected GA planet 
%         N.B (it's greater than the time needed to pass from v_inf_f to v_inf_i, 
%                     that is instead the time used to plot)
%   delta_deg      angle between v_inf_i and v_inf_f
%     hp[1]        pericenter altitude in [Km]

% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

ga_time=date2mjd2000(gaTime);

% DATA
mu_planet = astroConstants(10 + fb_Planet);                                % gravity constants of planet
G=astroConstants(1);
Mass_P=mu_planet/G;
R_planet = astroConstants(20 + fb_Planet);                                 % mean radius of planet
[kep_fb_Planet,mu_Sun] = uplanet(ga_time, fb_Planet);
Mass_sun=mu_Sun/G;
[r_planet, v_planet] = kep2car(kep_fb_Planet(1),kep_fb_Planet(2),kep_fb_Planet(3),kep_fb_Planet(4),kep_fb_Planet(5),kep_fb_Planet(6),mu_Sun);
r_SOI=norm(r_planet)*(Mass_P/Mass_sun)^(2/5);                              


% 2D HYPERBOLA DEFINITION
v_inf_i= V_minus-v_planet';
V_inf_i=norm(v_inf_i);
a_I=-mu_planet/V_inf_i^2;


v_inf_f= V_plus-v_planet';
V_inf_f=norm(v_inf_f);
a_f=-mu_planet/V_inf_f^2;

delta=acos(dot(v_inf_i, v_inf_f) /(V_inf_f*V_inf_i));                 % [rad]_Turning Angle
delta_deg=delta*180/pi;                                               % [deg]

rp_G=R_planet+300;
%options = optimset('TolX', 1e-13);
options = optimoptions('fsolve','FunctionTolerance',1e-13);

fun = @(rp) asin(mu_planet/(mu_planet+rp*V_inf_i^2))+asin(mu_planet/(mu_planet+rp*V_inf_f^2)) - delta;
%rp = fzero(fun, rp_G, options);                                      %[km]   rp_I==rp_f
rp = fsolve(fun, rp_G, options);


if rp<=rp_G
    rp=NaN;
end

hp=rp-R_planet;

e_I=1-rp/a_I;
e_f=1-rp/a_f;

DELTA_I=a_I^2*sqrt(e_I^2-1);
DELTA_f=a_f^2*sqrt(e_f^2-1);

if abs(DELTA_I+DELTA_f)-delta >= 0.1*pi/180
    disp 'error in Hyperbola_Plot'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYPERBOLAS ALWAYS LAYS IN SAME PLANE
vp_I=sqrt(2*mu_planet/rp+V_inf_i^2);                              
vp_f=sqrt(2*mu_planet/rp+V_inf_f^2);


% MISSION DEFINITION
n=cross(v_inf_i,v_inf_f)/norm(cross(v_inf_i,v_inf_f));

delta_I=2*asin(1/e_I);
[v_deltahalf]=v_rotated(v_inf_i, delta_I/2, n);
vp_hat=v_deltahalf/norm(v_deltahalf);

Vp_vec_I=vp_I*vp_hat;
Vp_vec_f=vp_f*vp_hat;
DV_cost_vec=Vp_vec_f-Vp_vec_I;

%DV_cost=abs(vp_f-vp_I);
%DV_cost_vec=DV_cost*vp_hat;

rp_hat=cross(vp_hat, n)/norm(cross(vp_hat, n));
rp_vec=rp*rp_hat;


time=1*24*3600;
N=600;
%%%% GET TIME TO GO FROM [-THETA_INF_I] TO [THETA=0]
TA_inf_I=acos(-1/e_I);
[TA_I, t] = KeplerSolverforHyperbola( e_I, a_I, mu_planet, 0, -TA_inf_I, time, N);
TA_inf_vec_I = TA_inf_I*ones(length(TA_I),1)';
sol = TA_inf_vec_I - TA_I;
index_I = find(sol<1e-10);
tt_I=t(index_I(1));
time_vec_I=linspace(0, tt_I, N);                                           %time from [-THETA_INF] to [+THETA_INF]


%%%% GET TIME TO GO FROM [THETA=0] TO [+THETA_INF_F] 
TA_inf_f=acos(-1/e_f);
[TA_f, t] = KeplerSolverforHyperbola( e_f, a_f, mu_planet, 0, -TA_inf_f, time, N);
TA_inf_vec_f = TA_inf_f*ones(length(TA_f),1)';
sol = TA_inf_vec_f - TA_f;
index_f = find(sol<1e-10);
tt_f=t(index_f(1));
time_vec_f=linspace(0, tt_f, N);                                           %time from [-THETA_INF] to [+THETA_INF]
  


initial_state_I=[rp_vec,  -Vp_vec_I];
initial_state_f=[rp_vec,  Vp_vec_f];

%[~, r_I, V_I_vec]=TPB_OrbitCompl(initial_state_I, mu_planet, time_vec_I/2);
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-13);

[~,state_I] = ode45(@(t,s) ode_2bp(t,s,mu_planet), time_vec_I/2, initial_state_I, ode_options);
r_I = [state_I(:,1), state_I(:,2),state_I(:,3)];
V_I_vec = [state_I(:,4);state_I(:,5);state_I(:,6)];

%[~, r_f, V_f_vec]=TPB_OrbitCompl(initial_state_f, mu_planet, time_vec_f/2);
[~,state_f] = ode45(@(t,s) ode_2bp(t,s,mu_planet), time_vec_f/2, initial_state_f, ode_options);
r_f = [state_f(:,1), state_f(:,2), state_f(:,3)];
V_f_vec = [state_f(:,4), state_f(:,5), state_f(:,6)];

Time_in_SOI=time_vec_I(N)/2 + time_vec_f(N)/2 + (r_SOI - norm(r_I(N,:)))/norm(V_I_vec(N,:)) + (r_SOI - norm(r_f(N,:)))/norm(V_f_vec(N,:));

%figure()
Venus_3D
%scatter3(0,0,0, 500, 'r', 'filled')
hold on
p1=plot3( r_I(:,1), r_I(:,2), r_I(:,3), 'b', 'linewidth', 1);
p2=scatter3 (r_I(1,1), r_I(1,2), r_I(1,3), 'k', 'linewidth', 0.1);
p3=plot3( r_f(:,1), r_f(:,2), r_f(:,3), 'r', 'linewidth', 1);
%p4=arrow3d([rp_vec(1), rp_vec(1)+1e5*DV_cost_vec(1)], [rp_vec(2), rp_vec(2)+1e5*DV_cost_vec(2)], [rp_vec(3), rp_vec(3)+1e5*DV_cost_vec(3)], 0.9, 0.5e3, 2*0.5e3, 'k');
p5=arrow3d( [0, 0.1e4*v_planet(1)], [0, 0.1e4*v_planet(2)], [0, 0.1e4*v_planet(3)], 0.9, 1e3, 2*1e3, 'g');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
legend ( [p1 p2 p3 p5], 'Incoming hyperbola', 'Pericenter', 'Outcoming hyperbola', 'Venus velocity direction')
title('Flyby in Venus-centred frame parallel to HECI')
grid on
end