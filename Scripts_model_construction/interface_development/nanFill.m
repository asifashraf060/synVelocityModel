function [nanFree_y] = nanFill(x, y)
%% This function removes nan values of 1-D y array based on nearest non-nan values
% INPUT:
%   x = 1-D array of numbers
%   y = 1-D array of numbers
%
% OUTPUT:
%   nanFree_y = y array filled with non-nan values in nan spaces

if length(x) == length(y)
    error(' Length of x and y must be equal')
end

% find the indices with nan values
ind1 = find(isnan(y) == 1);
% loop through the indices to make groups of nan indices
% this is done to linearly interpolate two end point values of a nan cluster
ds = [];
for i = 1:length(ind1)
    if i+1<=length(ind1)
    d       = abs(ind1(i)-ind1(i+1)); % to search for where the difference is more than one
    ds(:,i) = d;
    end
end
i_diff = find(ds>1); % indices where the difference is >1
                        % so those can be used to make groups
% loop to properly identify indices where 'ind1' jumps
b = [];                     % hint: Its both values around i_diff elements
for i = 1:length(i_diff)
    a = [i_diff(i) i_diff(i)+1];
    b(i,:) = a;
end
b = sort(b(:))';

%% FILL NAN INDICES WITH VALUES
if isempty(b) % if there is only 1 group of nans 
    value_fill = linspace(y(min(ind1)-1), y(max(ind1)+1), length(ind1));
    y(ind1)    = value_fill;
else          % if there are multiple groups of nans
    bb = ind1(b);
    
    cc = horzcat(ind1(1),bb, ind1(end)); % forcefully start and end
                                        % from 1st and last value of 'ind1'
    
    %loop through each groups of nan values and fill the y
    for i = 1:length(cc)-1
        
        %% SPECIAL CASE--BEGINNING OF ARRAY
        if cc(i) == 1 % it is special case, because there is prior element
            ll = y(cc(i+1)+1);
            value_fill = zeros(length(cc(i):1:cc(i+1)))+ll;
            value_fill = value_fill(1,:);
            y(cc(i):cc(i+1)) = value_fill;
        else
            %% NORMAL CASE--REST OF THE ARRAY
            kk = y(cc(i-1)); ll = y(cc(i+1)+1);
        end
        

        
    end
   
    end
end



end