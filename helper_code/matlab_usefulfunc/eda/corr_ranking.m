function [cc, pval, idx] = corr_ranking(datatable, column_num, corrtype)
%%
% compute correlation coefficient of one variable and all the others
%
% INPUT:
% datatable ... matrix of data (observations x variables)
% column_num ... the column number (variable of interest)
% corrtype ... 'Pearson' or 'Spearman'
%
% OUTPUT:
% cc ... vector of correlation coefficient
% pval ... corresponding p-values
% idx ... corresponding index of original columns
%
% EXAMPLE: [cc, pval, idx] = corr_ranking(rand(7, 30), 1, 'Pearson');
%

% inputs
if nargin < 2; column_num = 1; end
if nargin < 3; corrtype = 'Pearson'; end

% replace nan with median
datatable = nan_remove_pair(datatable, [], 'median');

% compute correlation coefficients
lenv = size(datatable, 2);
cc = nan(1, lenv);
pval = nan(1, lenv);
for c = 1:lenv
    [cc(c), pval(c)] = corr(datatable(:, column_num), datatable(:, c), 'type', corrtype);
end

% sort
[cc, idx] = sort(cc, 'descend');
pval = pval(idx);

% visualize
for i = 1:lenv - 1
    if pval(i+1) < 0.05
        col = 'r';
    else
        col = 'k';
    end
    bar(i, cc(i+1), 'FaceColor', col, 'EdgeColor', 'w')
    hold on;
end
xlabel('column index')
ylabel([corrtype ' correlation coefficient'])
set(gca, 'XTick', 1:(lenv-1), 'XTickLabel', idx(2:end))
set(gca, 'box', 'off', 'TickDir', 'out')