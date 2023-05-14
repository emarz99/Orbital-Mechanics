function [Gr_lon, Gr_midnight] = Get_Greenwich_longitude(date)
%--------------------------------------------------------------------------
%   Calculates the longitude of Greenwich meridian in a given date for year
%   2022.
%--------------------------------------------------------------------------
%   Form:
%   Gr_lon = Get_Greenwich_longitude(date)
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   date          (1x6)  Date vector in the same format as the input for date2mjd2000.    
%
%   -------
%   Outputs
%   -------
%   Gr_lon        [1]   Longitude of Greenwich (rad).
%
%--------------------------------------------------------------------------
% Programmed by: Jaime Fernández Diz
%
% Date:                  31/12/2022
% Revision:              
% Tested by:
%--------------------------------------------------------------------------

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

Gr_0 = 1.785385988925347;     %Force 0 at spring equinox
Gr_year = 2*pi*(sum(days_in_month(1:date(2)-1)) + date(3) - 1)/(366);
Gr_midnight = Gr_0 + Gr_year;
Gr_lon = Gr_midnight + 2*pi*(date(4)*3600 + date(5)*60 + date(6))/(24*3600);

Gr_lon = mod(Gr_lon,2*pi);

end