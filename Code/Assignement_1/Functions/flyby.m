function [DV_cost, DV_real, hp, delta] = flyby ( V_minus, V_plus, fb_Planet, fb_Date)

% flyby.m - Returns the cost of a GravityAssist/Fly-by maneuver 
%
% PROTOTYPE:
%   [DV_cost] = flyby ( V_minus, V_plus, fb_Planet, fb_Date)
%
% DESCRIPTION:
%   Returns a scalar value in [Km/s] which represents the norm of the 
%   DeltaV vector needed at pericenter  in order to pass from the 
%   Incoming to the Outcaming arc of hyperbola
%   N.B: in case DV_cost=0 ----> NOT POWERED GravityAssit
%
%
% INPUT:
%   V_minus[3]      Velocity vector [Km/s] containing the coordinates of velocity
%                   at the final point of the first Interplanetary-Leg
%   
%   V_plus[3]       Velocity vector [Km/s] containing the coordinates of velocity
%                   at the initial point of the second Interplanetary-Leg
%
%   fb_planet[1]    Integer number identifying the celestial body (< 11)
%                             1:   Mercury
%                             2:   Venus
%                             3:   Earth
%                             4:   Mars
%                             5:   Jupiter
%                             6:   Saturn
%                             7:   Uranus
%                             8:   Neptune
%                             9:   Pluto
%                             10:  Sun
%
%   fb_Date[1]     Time, modified Julian day since 01/01/2000, 12:00 noon
%                  (MJD2000 = MJD-51544.5)
%
%
% OUTPUT:
%  DV_cost[1]     Powered Fly-By Cost: DeltaV vector given at the
%                 hyperbola pericenter.
% 
%  DV_real[1]     Powered Fly-By Efficiency: norm of DeltaV between V_minus
%                 and V_plus
%
%   rp[1]         Perigee radius of Parabola [km]
%
%  delta[1]       Angle between V_minus and V_plus [rad]

% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

mu_planet = astroConstants(10 + fb_Planet);                                % gravity constants of planet
R_planet = astroConstants(20 + fb_Planet);                                 % mean radius of planet


V_planet = zeros(1,3);                                                     %[km/s]_orbital velocity in heliocentric frame
[kep_fb_Planet,mu_Sun] = uplanet(fb_Date, fb_Planet);
[~,V_planet(1,:)] = kep2car(kep_fb_Planet(1),kep_fb_Planet(2),kep_fb_Planet(3),kep_fb_Planet(4),kep_fb_Planet(5),kep_fb_Planet(6),mu_Sun);



% Compute the velocities relative to the planet before and after the fly-by
v_inf_i= V_minus- V_planet;
V_inf_i=norm(v_inf_i);



v_inf_f= V_plus-V_planet;
V_inf_f=norm(v_inf_f);


delta=acos(dot(v_inf_i, v_inf_f) /(V_inf_f*V_inf_i));                 % [rad]_Turning Angle

if fb_Planet == 1
rp_G=R_planet+200;
end
if fb_Planet == 2
rp_G=R_planet+300;
end


%options = optimset('TolX', 1e-13);
options = optimoptions('fsolve','FunctionTolerance',1e-13);

fun = @(rp) asin(mu_planet/(mu_planet+rp*V_inf_i^2))+asin(mu_planet/(mu_planet+rp*V_inf_f^2)) - delta;

%rp = fzero(fun, rp_G, options);                                   %[km]   rp_I==rp_f
rp = fsolve(fun, rp_G, options);


if rp < rp_G
    DV_cost=NaN;
    DV_real=NaN;
    hp=NaN;
    return
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYPERBOLAS ALWAYS LAYS IN SAME PLANE
DV_real=norm(v_inf_f-v_inf_i);
vp_I=sqrt(2*mu_planet/rp+V_inf_i^2);                              
vp_f=sqrt(2*mu_planet/rp+V_inf_f^2);
DV_cost=abs(vp_f-vp_I);
hp=rp-R_planet;
end

end