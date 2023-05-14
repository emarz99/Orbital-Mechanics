function PorkChop(earliest_dep_date, latest_arr_date, DV, N, toll, repT)

% flyby.m - Returns the Porckchop plot of the considered transfer
%
% PROTOTYPE:
%  PorkChop(earliest_dep_date, latest_arr_date, DV, N, toll, repT)
%
%
%
% INPUT:
%   earliest_dep_date [6]           Date in the Gregorian calendar, as a 6-element vector
%                                         [year, month, day, hour, minute, second]
%
%    latest_arr_date [6]            Date in the Gregorian calendar, as a 6-element vector
%                                         [year, month, day, hour, minute, second]
%
%         DV[NxN]                   Matrix of DV coming from Lambertsolver.m
%
%           N                       Numeber of points considered, needs to
%                                   be consistent with DV dimensions
%
%          toll                     tollerance around the minimum points
%                                   highlined on the plot, useful to 
%                                   graphically see repetitions patterns
%
% optionalINPUT: 
%          repT                     is the period of repetition graphically
%                                   computed from the porckchop, adding
%                                   this in the plot is highlibned a box
%                                   showing the zoomed repetition pattern,
%                                   it can also be used to zoom a certain
%                                   time window of length repT       

initialDepMJD2000 = date2mjd2000(earliest_dep_date);
latestArrMJD2000  = date2mjd2000(latest_arr_date);

MJD = linspace(initialDepMJD2000, latestArrMJD2000, N);

minDV = min(DV,[],'all');
[idx1, idx2] = find(DV == minDV);

tplot = NaN(N,1);
for i = 1:N
    dep_date=mjd20002date(MJD(i));
    tplot(i)= datenum(dep_date);   
end

tdep_minDV = tplot(idx1);
tarr_minDV = tplot(idx2);


all_minDV = NaN(N,1);
tdep_ALLminDV = NaN(N,1);
tarr_ALLminDV = NaN(N,1);

n = 1;

for j = 1:N
    for k = 1:N
        if abs(DV(j,k) - minDV) < toll
            
            all_minDV(n) = DV(j,k);
            tdep_ALLminDV(n) = tplot(j);
            tarr_ALLminDV(n) = tplot(k);
            n = n+1;
        end
    end
end

figure()
[C,h]=contour(tplot, tplot, DV',(1:1:8*minDV));
xlim([tplot(1) tplot(end)])
datetick('y', 'yyyy-mm-dd', 'keeplimits')
datetick('x', 'yyyy-mm-dd', 'keeplimits', 'keepticks')
xtickangle(45)
ytickangle(45)
grid on
hold on
caxis([minDV 3*minDV])
l = colorbar ;
ylabel(l, 'DV [km/s]');
clabel(C,h,(1:1:2*minDV)) 
xlabel('Departure date')
ylabel('Arrival date')
title('Porkchop plot for DV minimization')
scatter3(tdep_ALLminDV, tarr_ALLminDV, all_minDV,20,'filled','b')
scatter3(tdep_minDV, tarr_minDV, minDV,25,'filled','r')

if nargin == 6
    WindowINdays=latestArrMJD2000-initialDepMJD2000;
    stepINdays=WindowINdays/N;
    II=round((repT*365)/stepINdays);
    idxTransf=idx2-idx1;
    i=1;

    axes('position',[.20 .60 .25 .25])
    box on
    index_dep= (tplot > tplot(i)) & (tplot < tplot(i+II));
    index_arr= (tplot > tplot(i+idxTransf)) & (tplot < tplot(i+II+idxTransf));
    
    DVnew=DV(index_dep,index_arr);


    [C,h]=contour(tplot(index_dep), tplot(index_arr), DVnew',(1:1:8*minDV));
    %xlim([tplot(i) tplot(index_dep)])
    datetick('y', 'yyyy-mm-dd', 'keeplimits')
    datetick('x', 'yyyy-mm-dd', 'keeplimits', 'keepticks')
    xtickangle(45)
    ytickangle(45)
    grid on
end

end