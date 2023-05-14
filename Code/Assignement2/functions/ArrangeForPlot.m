function [x,y] = ArrangeForPlot(x,y)

% Inserts NaN in a set of longitudes and latitudes or right ascensions or
% declinations in order to avoid horizontal lines in plots. 
%
% PROTOTYPE:
% [x, y] = ArrangeForPlot(x, y)
%
% INPUT:
% x(1,:) Longitude or right ascension [-][rad]
% y(1,:) Latitude or declination [-][rad]
% OUTPUT:
% x(1,:) Longitude or right ascension [-][rad]
% y(1,:) Latitude or declination [-][rad]
%
% CONTRIBUTORS:
% Jaime Fernández Diz
%
% VERSIONS
% 2022 12 06: First version

for i = length(x):-1:2
    if x(i) < 0 && x(i-1)>= 0.8*pi
        aux = x(i:end);
        x(i+1:end+1) = aux;
        x(i) = NaN;
        aux = y(i:end);
        y(i+1:end+1) = aux;
        y(i) = NaN;
    elseif x(i) > 0 && x(i-1)<= -0.8*pi
        aux = x(i:end);
        x(i+1:end+1) = aux;
        x(i) = NaN;
        aux = y(i:end);
        y(i+1:end+1) = aux;
        y(i) = NaN;
    end
end


end