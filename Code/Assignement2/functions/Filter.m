function filtered = Filter(time_serie, npoints)

%--------------------------------------------------------------------------
%   Applies movemean to a set of data in columns.
%--------------------------------------------------------------------------
%   Form:
%   secular = Filter(keplerian_elements, 230)
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   time_serie    (:,:)   Set of data in columns. [Any units]
%   npoints       [1]     Number of points to be averaged. [#]      
%
%   -------
%   Outputs
%   -------
%   filtered      (:,:)   Set of data in columns. [Any units]
%
%--------------------------------------------------------------------------
% Programmed by: Jaime Fernández Diz
%
% Date:                  24/12/2022
% Revision:              
% Tested by:
%--------------------------------------------------------------------------
filtered = zeros(size(time_serie));

for ii = 1:size(time_serie,2)
    filtered(:,ii) = movmean(time_serie(:,ii), npoints);
end


end