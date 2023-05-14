function [f_co] = lowPassFilter(vec)

% lowPassFilter - Select a proper cut-off frequency for a non-monothonic
% oscillations. The cut-off frequency is selected studying how many points
% are in an ascending or descending branch

% INPUT:
% vec       [1xN]   vector containing the function values at each instant
%                   considered

% OUTPUT:
% f_co      [1x1]   cut-off frequency [Hz]

% AUTHORS:
% Pasquariello Chiara
% Ferro Jacopo
% Giorgini Francesco
% Guidetti Tommaso


K = [];
k = 1; % initial position index in the vector

% To search if the branch is ascending or descending through the vector

while k < length(vec)

    if vec(k) <= vec(k+1)
% Acending branch
        n = 1;
    else
% Descending branch        
        n = 2;
    end

% To count how many points are in the ascending or descending branch

    switch n
        case 1
            i = 0; % initialize the index that counts the points in the ascending branch
            while vec(k) < vec(k+1)
                k = k+1; % position index in the vector
                i = i+1; % index that counts the points in the ascending branch
                if k == length(vec) % when the last vector element is reached the last 
                                    % branch is deleted from the computation
                    i = NaN;
                    break
                end
            end
            K = [K;i]; % To save the number of points of each branch
        case 2
            i = 0; % initialize the index that counts the points in the descending branch
            while vec(k) > vec(k+1)
                k = k+1; % position index in the vector
                i = i+1; % index that counts the points in the descending branch
                if k == length(vec) % when the last vector element is reached the last 
                                    % branch is deleted from the computation
                    i = NaN;
                    break
                end
            end
            K = [K;i]; % To save the number of points of each branch
    end

end

% To compute the cut-off frequency

f_co = 2*round(mean(K,'omitnan'));

end


