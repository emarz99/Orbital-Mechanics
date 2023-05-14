function [a3x, a3y, a3z] = thirdbody(mu, x2, y2, z2, x3, y3, z3)
%--------------------------------------------------------------------------
%   Computes the acceleration induced by a body different from the main
%   one, taking on account the direct and indirect perturbation.
%   Input may be a vector. In that case, each of the terms correspond to a 
%   different body.
%--------------------------------------------------------------------------
%   Form:
%   [a3x, a3y, a3z] = thirdbody(mu, x, y, z);
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   mu            (:,1)   Gravitational parameter of each body.
%   x2, y2, z2    (1)     Position of the body whose orbit is being propagated
%   x3, y3, z3    (:,1)   Position of the perturbator body(es) on the chosen frame.             
%
%   -------
%   Outputs
%   -------
%   a3x, a3y, a3z    (1)   Acceleration on the chosen frame.
%
%--------------------------------------------------------------------------
% Programmed by: Jaime Fernández Diz
%
% Date:                  08/03/2020
% Revision:              
% Tested by:
%--------------------------------------------------------------------------

r23 = zeros(max(size(mu,1), size(mu,2)),1);
r3 = r23;
a3xi = r23;
a3yi = r23;
a3zi = r23;

for i = 1:max(size(mu,1), size(mu,2))    
    r23(i) = sqrt( (x2-x3(i))^2 + (y2 - y3(i))^2 + (z2 - z3(i))^2 );
    r3(i) = sqrt( x3(i)^2 + y3(i)^2 + z3(i)^2 );
    a3xi(i) = -mu(i)*(-(x3(i)-x2)/(r23(i)^3) + x3(i)/r3(i)^3);
    a3yi(i) = -mu(i)*(-(y3(i)-y2)/(r23(i)^3) + y3(i)/r3(i)^3);
    a3zi(i) = -mu(i)*(-(z3(i)-z2)/(r23(i)^3) + z3(i)/r3(i)^3);
    
end

a3x = sum(a3xi);
a3y = sum(a3yi);
a3z = sum(a3zi);

