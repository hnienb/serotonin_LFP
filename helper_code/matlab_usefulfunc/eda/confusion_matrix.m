function [cm, pm] = confusion_matrix(datatable, corrtype, mapcolor, textoption)
% generate 
% confusion matrix (correlation coefficient between pairs of the
% given vectors)
%
% INPUT:
% datatable ... matrix of data (observations x variables)
% corrtype ... 'Pearson' or 'Spearman'
% mapcolor ... color scheme of confusion matrix ('parula' in default)
% textoption ... 0, no text on the colormap; 1, add p-val
%
% OUTPUT:
% cm ... confusion matrix
% pm ... matrix for corresponding p-values
%

% inputs
if nargin < 2; corrtype = 'Pearson'; end
if nargin < 3; mapcolor = 'jet'; end
if nargin < 3; textoption = 1; end

% replace nan with median
datatable = nan_remove_pair(datatable, [], 'median');

% generate confusion matrix
lenv = size(datatable, 2);
cm = nan(lenv, lenv);
pm = nan(lenv, lenv);
for r = 1:lenv
    for c = 1:lenv
        [cm(r,c), pm(r,c)] = corr(datatable(:, r), datatable(:, c), 'type', corrtype);
    end
end

% visualize
imagesc(1:lenv, 1:lenv, cm)
colormap(mapcolor)
colorbar('eastoutside')
if textoption
    for r = 1:lenv
        for c = 1:lenv
            if pm(r, c) < 0.05
                col = 'w';
            else
                col = 'k';
            end
            text(r-0.4, c, num2str(pval_inequality(pm(r, c))), 'color', col)
        end
    end
end
set(gca, 'box', 'off', 'TickDir', 'out')