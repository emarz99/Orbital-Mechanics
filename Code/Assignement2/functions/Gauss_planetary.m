function ds = Gauss_planetary(t, KEP_ELEMS, mu_planet, perturbations)

%--------------------------------------------------------------------------
%   Computes the keplerian state of a body through Gauss planetary
%   equations in TNH frame.
%--------------------------------------------------------------------------
%   Form:
%   ds = Gauss_planetary(t, KEP_ELEMS, mu_earth, @(t,s) a_moon(t,s,data))
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   t             [1]     Instant of time of the propagation. [T][s]
%   KEP_ELEMS     (1,6)   Keplerian state (SMA, ECC, INC, AOP, RAAN, TA). 
%                         [L, ~, ~, ~, ~, ~][Km,-, rad, rad, rad, rad].
%   mu_planet     [1]     Gravitational parameter of the central body. 
%                         [L^3/M T^2][km^3/(kg*s^2)].
%   perturbations [1]     Pointer to a function which returns the
%                         perturbation acceleration in cartesian 
%                         coordinates.
%   -------
%   Outputs
%   -------
%   ds            (:,6)   Derivative of the keplerian state (SMA, ECC, INC, 
%                         AOP, RAAN, TA). [L/T, ~, ~, ~, ~, ~][Km/s, 1/s, 
%                         rad/s, rad/s, rad/s, rad/s]
%
%--------------------------------------------------------------------------
% Programmed by: Jaime Fernández Diz
%
% Date:                  20/12/2022
% Revision:              
% Tested by:
%--------------------------------------------------------------------------
SMA = KEP_ELEMS(1);
ECC = KEP_ELEMS(2);
INC = KEP_ELEMS(3);
AOP = KEP_ELEMS(4);
RAAN = KEP_ELEMS(5);
TA = KEP_ELEMS(6);

[r_vec, v_vec] = kep2car(SMA,ECC,INC*180/pi,RAAN*180/pi,AOP*180/pi,TA*180/pi,mu_planet);
s_car = [r_vec, v_vec];

u_t = v_vec/norm(v_vec);
u_h = cross(r_vec, v_vec)/norm(cross(r_vec, v_vec));
u_n = cross(u_h, u_t);

a_p_xyz = perturbations(t, s_car);
a_p = [dot(a_p_xyz, u_t), dot(a_p_xyz, u_n), dot(a_p_xyz, u_h)];

p = SMA*(1 - ECC^2);
r = p/(1 + ECC*cos(TA));
vel = sqrt(2*mu_planet/r - mu_planet/SMA);
h = sqrt(p*mu_planet);

da = a_p(1)*2*SMA^2*vel/mu_planet;
de = (2*(ECC + cos(TA))*a_p(1) - r*sin(TA)*a_p(2)/SMA)/vel;
di = r*cos(TA + AOP)*a_p(3)/h;

dw = (2*sin(TA)*a_p(1) + (2*ECC + r*cos(TA)/SMA)*a_p(2))/(ECC*vel) - a_p(3)*r*sin(TA + AOP)*cos(INC)/(h*sin(INC));
dW = a_p(3)*r*sin(TA + AOP)/(h*sin(INC));
dtheta = h/r^2 - (2*sin(TA)*a_p(1) + (2*ECC + r*cos(TA)/SMA)*a_p(2))/(ECC*vel);

ds = [da, de, di, dw, dW, dtheta]';

end