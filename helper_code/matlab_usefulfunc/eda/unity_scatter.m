function unity_scatter(x, y, label, visparams)
%%
% generate a scatter plot with the unity line
% INPUT:
% x, y ... vector with the same length
% label ... interger vector with the same length, indicating a group 
% visparams ... struct containing parameters for appearance
%                 EXAMPLE: 
%                 visparams.xynames = {'variable A', 'variable B'}
%                 visparams.groupnames = {'animal M', 'animal K'};
%                 visparams.marker = {'o', 's'};
%                 visparams.markersize = 50;
%                 visparams.markerfacealpha = 0.4;
%                 visparams.markeredgealpha = 0.6;
%                 visparams.colors = [0 0.5 0; 1 0 0];
%                 visparams.range = [0 1];
%                 visparams.fontsize = 8;
%
% EXAMPLE: unity_scatter(randn(30,1), randn(30,1), randi(2, 30, 1)-1)
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% check legnth of x, y
if size(x, 1)==1; x = x'; end
if size(y, 1)==1; y = y'; end
if nargin < 2 && size(x, 2)==2; y = x(:, 2); x = x(:, 1); end
if nargin < 3; label = ones(size(x)); end
if nargin < 4
    visparams.xynames = {'', ''};
end
if length(x)~=length(y)
    error('vector x and y must be the same-length')
end
if length(label)~=length(x)
    error('label length must be the same as the data')
end

% parameters default
if ~isfield(visparams, 'xynames')
    visparams.xynames = {'', ''};
end
if ~isfield(visparams, 'groupnames')
    visparams.groupnames = [];
end
if ~isfield(visparams, 'marker')
    visparams.marker = [];
end
if ~isfield(visparams, 'markersize')
    visparams.markersize = 30;
end
if ~isfield(visparams, 'markerfacealpha')
    visparams.markerfacealpha = 0.4;
end
if ~isfield(visparams, 'markeredgealpha')
    visparams.markeredgealpha = 0.6;
end
if ~isfield(visparams, 'colors')
    visparams.colors = [];
end
if ~isfield(visparams, 'range')
    visparams.range = [];
end
if ~isfield(visparams, 'fontsize')
    visparams.fontsize = 8;
end

% range
all = [x(:); y(:)];
dist = max(all) - min(all);
if isempty(visparams.range)
    visparams.range = [min(all) - 0.1*dist, max(all) + 0.1*dist];
end

% unity
plot(visparams.range, [0, 0], '-','color',0.4*ones(1,3))
hold on;
plot([0, 0], visparams.range, '-','color',0.4*ones(1,3))
hold on;
plot(visparams.range, visparams.range, '-','color',0.4*ones(1,3))

% scatter
groups = unique(label);
n_group = length(groups);
if n_group==1
    stats = cell(1, 1);
    if isempty(visparams.colors)
        visparams.colors = zeros(1, 3);
    end
    if isempty(visparams.marker)
        visparams.marker = {'o'};
    end
    visparams.groupnames = '';
else
    stats = cell(1, n_group + 1);
    if isempty(visparams.colors)
        visparams.colors = lines(n_group);
    end
    if size(visparams.colors, 1) == n_group
        visparams.colors = [visparams.colors; 0 0 0];
    end
    if isempty(visparams.groupnames)
        visparams.groupnames = cell(1, n_group);
        for n = 1:n_group
            visparams.groupnames{n} = ['group ' num2str(n)];
        end
    end
    if length(visparams.marker) < n_group
        if isempty(visparams.marker)
            marker = 'o';
        else
            marker = visparams.marker;
        end
        visparams.marker = cell(1, n_group);
        for n = 1:n_group
            visparams.marker{n} = marker;
        end
    end
end

% plot
for n = 1:n_group
    hold on;
    xd = x(label==groups(n)); yd = y(label==groups(n));
    [xd, yd] = nan_remove_pair(xd, yd);
    scatter(xd, yd, visparams.markersize, 'filled', 'marker', visparams.marker{n}, ...
        'markerfacecolor', visparams.colors(n,:), ...
        'markerfacealpha', visparams.markerfacealpha, ...
        'markeredgecolor', 'w', 'markeredgealpha', visparams.markeredgealpha)
    stats{n} = pair_tests([xd, yd]);    
end

% axis
axis([visparams.range visparams.range])
set(gca, 'box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [visparams.range(1) visparams.range(end)])
set(gca, 'YTick', [visparams.range(1) visparams.range(end)])

% stats across groups
if n_group > 1
   hold on;
   xd = x; yd = y;
   [xd, yd] = nan_remove_pair(xd, yd);
   stats{n+1} = pair_tests([xd, yd]);
end

% display stats (non-parametric paired-test)
lens = length(stats);
for n = 1:lens
    if n==lens
        g_lab = 'all';
    else
        g_lab = visparams.groupnames{n};
    end
    if stats{n}.pair(1).signrank.p < 0.05
        text(visparams.range(1)+0.1*(visparams.range(2)-visparams.range(1)), ...
            visparams.range(1)+(0.975 - 0.075*(n-1))*(visparams.range(2)-visparams.range(1)), ...
            [g_lab ': n=' num2str(stats{n}.pair(1).n) ', p<' num2str(pval_inequality(stats{n}.pair(1).signrank.p))], ...
            'fontsize', visparams.fontsize, 'color', visparams.colors(n,:))
    else
        text(visparams.range(1)+0.1*(visparams.range(2)-visparams.range(1)), ...
            visparams.range(1)+(0.975 - 0.075*(n-1))*(visparams.range(2)-visparams.range(1)), ...
            [g_lab ': n=' num2str(stats{n}.pair(1).n)...
            ', p=' num2str(pval_inequality(stats{n}.pair(1).signrank.p))], ...
            'fontsize', visparams.fontsize, 'color', visparams.colors(n,:))
    end
end        
% display stats (non-parametric across group test)
if n_group > 1
    C = nchoosek(1:n_group, 2);
    sz = size(C);
    col = pink(sz(1));
    for i = 1:sz(1)
        xd = x(label==groups(C(i, 1))); yd = y(label==groups(C(i, 2)));
        xd = nan_remove_pair(xd, []);
        yd = nan_remove_pair(yd, []);
        pval = ranksum(xd, yd)*sz(1); % bonferroni correction
        if pval < 0.05
            xr = (visparams.range(1)+(0.1 - i*0.03)*(visparams.range(2)-visparams.range(1))).*ones(1, 2);
            yr = [visparams.range(1)+(0.975 - 0.075*(C(i, 1)-1))*(visparams.range(2)-visparams.range(1)), ...
                visparams.range(1)+(0.975 - 0.075*(C(i, 2)-1))*(visparams.range(2)-visparams.range(1))];
            hold on;
            plot(xr, yr, '-', 'color', col(i, :), 'linewidth', 0.5)
            hold on;
            plot(xr(1), mean(yr), '*', 'color', col(i, :))
        end
    end
end
xlabel(visparams.xynames{1}, 'fontsize', visparams.fontsize)
xlabel(visparams.xynames{2}, 'fontsize', visparams.fontsize)
% axis square

