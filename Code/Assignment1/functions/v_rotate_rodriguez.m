function [v_rotated] = v_rotate_rodriguez(v,plane,angle)
% This function performs Rodriguez' Formula to rotate
% vector [v] around vector [plane vector] certain angle following right
% hand rule.
% 
% PROTOTYPE
%     v_new = v_rotate_rodriguez(v,plane,angle)
% 
% INPUT:
%     v[3]          Vector to be rotated [L]
%     plane[3]      Normal plane vector [L]
%     angle[1]      Angle to be rotated [rad]
%
% OUTPUT:
%     v_rotated[3]  Rotated vector [L]  
%
% CONTRIBUTORS:
%     Carlos Albiñana
%
% VERSIONS:
%     2022-12-9:    First Version


v_rotated=v*cos(angle)+(cross(plane,v)*sin(angle))+plane*dot(plane,v)*(1-cos(angle));


end