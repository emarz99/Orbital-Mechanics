function [TA,t] = KeplerSolverforHyperbola( e, a, mu, t0, TA0, time, N )

% KeplerSolverforHyperbola.m - TA (True Anomaly) evolution and the t (time)     
%                              associated
%                            
% PROTOTYPE:
%   [TA,t] = KeplerSolverforHyperbola( e, a, mu, t0, TA0, time, N )
%
% DESCRIPTION:
%   Returns a vector of TA and a same dimension t vector associated,
%   hyperbolic function, Mean and Eccentric anomaly are used to compute TA
%   of hyperbola at each time istant t
%
%
% INPUT:
%   e[1]       eccentricity [-] of hyperbola     (>1)
%   a[1]       semimajor axis [km] of hyperbola  (< 0)
%   mu[1]      planetetary gravity constant [Km^3/s^2] of the planet acting as focus
%              of hyperbola
%  t0[1]       initial value of time vector (t) considered
%  TA0[1]      TA value at t0
%  time[1]     final value of time vector (t) considered
%   N[1]       lenght of output vectors 
%
%
% OUTPUT:
%  TA[1xN]     True anomaly vector in the time window specified in input
%  t[1xN]      Time vector composed by istants in which TA is calculated

% AUTORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso

t = linspace (t0, time, N);
F0 = 2*atan(sqrt((e-1)/(1+e))*tanh(TA0/2));
M0 = F0-e*sinh(F0);
TA=zeros(1,length(t));
for i=1:length(t)
    M = M0+sqrt(mu/((-a)^3))*(t(i)-t0);
    Fg = M+(e*sinh(M)/(1-sinh(M+e)+sinh(M)));

    k=isnan(Fg);
    if k==1
        Fg=0;
    end


    fun = @(F) M+F-e*sinh(F);
   
    options = optimoptions('fsolve','FunctionTolerance',1e-13);
    F = fsolve(fun,Fg,options);

%     options = optimset('TolX',1e-16);
%     F = fzero(fun,Fg,options);
    

    f = 2*atan(sqrt((1+e)/(e-1))*tanh(F/2));
    if f<0
        f = f+2*pi;
    end
    TA(i) = f;
end
end