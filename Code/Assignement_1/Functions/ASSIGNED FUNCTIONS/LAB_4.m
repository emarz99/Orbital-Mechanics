
%% DATA

mu_E=astroConstants(13);
mu_S=astroConstants(4);
AU=astroConstants(2);

DELTA=9200;                       %[km]
r_E=[1, 0, 0]*AU;                 %[km]
v_inf_i=[ 15.1, 0, 0 ];           %[km/s]

%2D HYPERBOLA
v_inf=norm(v_inf_i);               % (v_inf-=v_inf+) in magnitude
a=-mu_E/(v_inf^2);

delta=2*atan(-a/DELTA);            % [rad]_Turning Angle 
e=1/sin(delta/2);                  
rp= a*(1-e);

v_e=sqrt(mu_S/norm(r_E));
V_e=[0, v_e, 0];

n1=[0, 0, 1];    % Behind planet
n2=[0, 0, -1];   % In Front of Planet
n3=[0, -1, 0];   % Under the planet
n4=[0, -1, 0];   % Up the planet

