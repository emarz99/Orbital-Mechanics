function rho = expAtmModel(h)

%expAtmModel: returns the value of air density computed with the exponential
%             atmospheric model. 
%             h0      Reference altitude [km]
%             rho0    Reference density [kg/m^3]
%             H       Scale Factor [km]

% INPUT:
% h       [1x1] Altitude of the satellite [km]

% OUTPUT:
% rho     [1x1] Air density [kg/m^3]

% AUTORS:
% Ferro Jacopo
% Giorgini Francesco
% Guidetto Tommaso
% Pasquariello Chiara

if h >=0 && h < 25
    h0 = 0;
    rho0 = 1.225;
    H = 7.249;
elseif h >=25 && h < 30
    h0 = 25;
    rho0 = 3.899*1e-2;
    H = 6.349;
elseif h >=30 && h < 40
    h0 = 30;
    rho0 = 1.774*1e-2;
    H = 6.682;
elseif h >=40 && h < 50
    h0 = 40;
    rho0 = 3.972*1e-3;
    H = 7.554;
elseif h >=50 && h < 60
    h0 = 50;
    rho0 = 1.057*1e-3;
    H = 8.382;
elseif h >=60 && h < 70
    h0 = 60;
    rho0 = 3.206*1e-4;
    H = 7.714;
elseif h >=70 && h < 80
    h0 = 70;
    rho0 = 8.770*1e-5;
    H = 6.5799;
elseif h >= 80 && h < 90
    h0 = 80;
    rho0 = 1.905*1e-5;
    H = 5.799;
elseif h >= 90 && h < 100
    h0 = 90;
    rho0 = 3.396*1e-6;
    H = 5.382;
elseif h >= 100 && h < 110
    h0 = 100;
    rho0 = 5.297*1e-7;
    H = 5.877;
elseif h >= 110 && h < 120
    h0 = 110;
    rho0 = 9.661*1e-8;
    H = 7.263;
elseif h >= 120 && h < 130
    h0 = 120;
    rho0 = 2.438*1e-8;
    H = 9.473;
elseif h >= 130 && h < 140
    h0 = 130;
    rho0 = 8.484*1e-9;
    H = 12.636;
elseif h >= 140 && h < 150
    h0 = 140;
    rho0 = 3.845*1e-9;
    H = 16.149;
elseif h >= 150 && h < 180
    h0 = 150;
    rho0 = 2.070*1e-9;
    H = 22.523;
elseif h >= 180 && h < 200
    h0 = 180;
    rho0 = 5.464*1e-10;
    H = 29.740;
elseif h >= 200 && h < 250
    h0 = 200;
    rho0 = 2.789*1e-10;
    H = 37.105;
elseif h >= 250 && h < 300
    h0 = 250;
    rho0 = 7.248*1e-11;
    H = 45.546;
elseif h >= 300 && h < 350
    h0 = 300;
    rho0 = 2.418*1e-11;
    H = 53.628;
elseif h >= 350 && h < 400
    h0 = 350;
    rho0 = 9.158*1e-12;
    H = 53.298;
elseif h >= 400 && h < 450
    h0 = 400;
    rho0 = 3.725*1e-12;
    H = 58.515;
elseif h >= 450 && h < 500
    h0 = 450;
    rho0 = 1.585*1e-12;
    H = 60.828;
elseif h >= 500 && h < 600
    h0 = 500;
    rho0 = 6.967*1e-13;
    H = 63.822;
elseif h >= 600 && h < 700
    h0 = 600;
    rho0 = 1.454*1e-13;
    H = 71.835;
elseif h >= 700 && h < 800
    h0 = 700;
    rho0 = 3.614*1e-14;
    H = 88.667;
elseif h >= 800 && h < 900
    h0 = 800;
    rho0 = 1.170*1e-14;
    H = 124.64;
elseif h >= 900 && h < 1000
    h0 = 900;
    rho0 = 5.245*1e-15;
    H = 181.05;
elseif h >= 1000
    h0 = 1000;
    rho0 = 3.019*1e-15;
    H = 268.00;
else 
    disp('Value of h not assigned')
end

rho = rho0*exp(-(h-h0)/H);

end