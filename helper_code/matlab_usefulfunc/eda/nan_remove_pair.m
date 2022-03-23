function [x, y] = nan_remove_pair(x, y, replace)
%%
% remove nans from paired-data
%
% INPUT: x, y ... vectors with the same length
%        replace ... 'none' (remove nans)
%                    'mean' (replace nans with mean)
%                    'median' (replace nans with median)
%
%        when only x (matrix) is given, the function performs nan_removal for all the rows
%
% EXAMPLE: [x, y] = nan_remove_pair([nan; 3; 4], [2; nan; 5]);
%          x = nan_remove_pair([nan; 3; 4, 2; nan; 5, 1; 2; 3], [], 'none');
%

if nargin < 3; replace = 'none'; end

szx = size(x, 2);
if szx == 1
    X = [x, y];
else
    X = x;
end
switch replace
    case 'none'
        nans = any(isnan(X),2);
        X(nans, :) = []; 
    case 'mean'
        for i = 1:size(X, 2)
            X(isnan(X(:, i)), i) = nanmean(X(:, i));
        end
    case 'median'
        for i = 1:size(X, 2)
            X(isnan(X(:, i)), i) = nanmedian(X(:, i));
        end 
end
if szx == 1
    x = X(:, szx);
    y = X(:, szx+1:end);
else
    x = X;
end