function [v_inf_f]=v_rotated(v_inf_i, delta, u)

% v_rotated.m - computes the rotation of vectors
%
% PROTOTYPE
% [v_inf_f]=v_rotated(v_inf_i, delta, u)
%
% DESCRIPTION: rotates the input velocity in a counterclocwise direction
%               around a versor U of an angle DELTA
%
% INPUT:
%    v_inf_i[3]     input  vector in 
%    delta[1]       angle in [rad]
%      u[3]         normal versor individuating the plane on which v_inf_i
%                   lays and then the rotation direction
%
%
% OUTPUT:
%    v_inf_f[3]     is the input vector rotated of delta so has
%                   the same magnitude and different direction

% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

v=v_inf_i;
v_rotated=v*cos(delta)+cross(u,v)*sin(delta)+u*(dot(u,v))*(1-cos(delta));
v_inf_f=v_rotated;
end