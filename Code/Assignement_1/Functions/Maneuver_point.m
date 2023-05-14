function[p1, p2, p3, p4, p5]=Maneuver_point(V1, V2, Origin, V_3)

% Maneuver_point.m - plot in 3Dspace of the Velocity Triangle of a generic
%                    maneuver point
%
% PROTOTYPE:
%    [p1, p2, p3, p4, p5]=Maneuver_point(V1, V2, Origin, V_fb_p)
%
% DESCRIPTION: the Velocity Triangle is composed of three vectors:
%                          -Velocity before maneuvre
%                          -Velocity after maneuver
%                          -Difference of the previous vectors
%
%
% INPUT:
%   V1[3]     Velocity vector [Km/s] Before the Maneuvre
%   V2[3]     Velocity Vector [Km/s] After the Maneuvre
%  Origin[3]  Initial point of V1 & V2 coordinates (Maneuvre point)
%                      (default: [0, 0, 0])
%
% optional INPUT:
% V_3[3]      Velocity vector at pericenter of hyperbola (it is the cost of
%             the maneuvre but it,s not the conjunction of Velocity vectors
%             before and after GA): plotted to compare.
%
%
%OUTPUT: 
%  p1         Origin point plot
%  p2         V1 arrow-plot
%  p3         V2 arrow-plot
%  p4         DV=V2-V1 arrow plot
%
%optional OUTPUT:
%  p5          V_fb_p arrow plot
%
%
%NOTE:         output are needed to define appropriate legend outside the
%              function (example under):
%     legend( [p1, p2, p3, p4, p5] , ' maneuver point', 'Velocity before maneuver', 'Velocity after Maneuver', 'Maneuver effect', 'Maneuver Impulse') 
%
% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso


head_frac=0.9;
radii=0.1;
radii2=2*radii;

x1=[Origin(1), Origin(1)+V1(1)];
y1=[Origin(2), Origin(2)+V1(2)];
z1=[Origin(3), Origin(3)+V1(3)];

x2=[Origin(1), Origin(1)+V2(1)];
y2=[Origin(2), Origin(2)+V2(2)];
z2=[Origin(3), Origin(3)+V2(3)];

DX=[Origin(1)+V1(1), Origin(1)+V2(1)];
DY=[Origin(2)+V1(2), Origin(2)+V2(2)];
DZ=[Origin(3)+V1(3), Origin(3)+V2(3)];


figure()
hold on
p1=scatter3( Origin(1), Origin(2), Origin(3), 50, 'filled', 'ko');
p2=arrow3d(x1,y1,z1, head_frac,radii,radii2, 'r');
p3=arrow3d(x2,y2,z2, head_frac,radii,radii2, 'b');
p4=arrow3d(DX, DY, DZ,head_frac,radii,radii2,'g');


if nargin==4
    x3=[Origin(1)+V1(1), Origin(1)+V1(1)+V_3(1)];
    y3=[Origin(2)+V1(2), Origin(2)+V1(2)+V_3(2)];
    z3=[Origin(3)+V1(3), Origin(3)+V1(3)+V_3(3)];
    p5=arrow3d(x3, y3, z3, head_frac,radii,radii2, 'k');
end

% lim=max([x1, y1, z1, x2, y2, z2, DX, DY, DZ]);
% xlim([-lim, lim]);
% ylim([-lim, lim]);
% zlim([-lim, lim]);

xlabel(' X [km/s]')
ylabel(' Y [km/s]')
zlabel(' Z [km/s]')
grid on
axis equal
end

% legend( [p1, p2, p3, p4, p5] , ' maneuver point', 'Velocity before maneuver', 'Velocity after Maneuver', 'Maneuver effect', 'Maneuver Impulse')
% title('')


