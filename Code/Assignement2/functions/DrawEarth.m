function [] = DrawEarth()

% Draws the Earth in the current active figure
%
% PROTOTYPE:
% DrawEarth()
%
% INPUT:
%
% OUTPUT:
%
% CONTRIBUTORS:
% Jaime Fernández Diz
%
% VERSIONS
% 2022 09 21: First version

axis equal
[Xt,Yt,Zt] = sphere(); %Creates an sphere wit 100 spherical blocks 
Xt = Xt * 6371; %You multiply the coordinates times the Earth radius
Yt = Yt * 6371; %You multiply the coordinates times the Earth radius
Zt = Zt * 6371; %You multiply the coordinates times the Earth radius
[atlante] = imread('earth_globe.JPG'); %Used to print the image in the surface
Globe = surf(Xt,Yt,Zt); %Creates the surface
hT = findobj(Globe,'Type','Surface'); %You specify that the object Globe is a surface
set(hT,'CData',atlante,'facecolor','texturemap','edgecolor','none'); %You add the image to the surface
end