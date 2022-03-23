function ls_scatter(x, y, label, permutation_option)
%%
% generate a scatter plot with a linear regression line
%
% INPUT:
% x, y ... vector with the same length
% label ... interger vector with the same length, indicating a group 
% permutation_option ... 0 or 1; perform a permutation test to see if the correlation
% coefficients between the two paired groups are significantly different
% (this option only works when there are two unique labels in the 'label')
%
% EXAMPLE: ls_scatter(randn(30,1), randn(30,1), randi(2, 30, 1)-1, 1)
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% check legnth of x, y
if size(x, 1)==1; x = x'; end
if size(y, 1)==1; y = y'; end
if nargin < 2 && size(x, 2)==2; y = x(:, 2); x = x(:, 1); end
if nargin < 3; label = ones(size(x)); end
if nargin < 4; permutation_option = 0; end
if length(x)~=length(y)
    error('vector x and y must be the same-length')
end
if length(label)~=length(x)
    error('label length must be the same as the data')
end

% range
dist = max(x(:)) - min(x(:));
xrange = [min(x(:)) - 0.1*dist, max(x(:)) + 0.1*dist];
xrange(1) = floor(10*xrange(1))/10;
xrange(2) = ceil(10*xrange(2))/10;
dist = max(y(:)) - min(y(:));
yrange = [min(y(:)) - 0.1*dist, max(y(:)) + 0.1*dist];
% yrange(1) = floor(10*yrange(1))/10;
% yrange(2) = ceil(10*yrange(2))/10;

% scatter
groups = unique(label);
n_group = length(groups);
map = lines(n_group);
stats = nan(n_group, 5);
for n = 1:n_group
    hold on;
    
    % nan removal
    xd = x(label==groups(n)); yd = y(label==groups(n));
    [xd, yd] = nan_remove_pair(xd, yd);
    stats(n, 5) = length(xd);
    
    % scatter
    scatter(xd, yd, 20, 'filled','marker','o', 'markerfacecolor', map(n,:), ...
        'markerfacealpha',0.4, 'markeredgecolor','w','markeredgealpha',0.8)
    
    % linear regression
    beta = glmfit(xd, yd);
    hold on;
    plot(xrange, glmval(beta, xrange, 'identity'), '-', 'color', map(n, :))
    
    % correlation coefficients
    [rr, pp] = corrcoef(xd, yd);
    stats(n, 1:2) = [rr(1,2), pp(1,2)];
%     [stats(n, 3), stats(n, 4)] = corr(xd, yd, 'type', 'Spearman');
    [stats(n, 3), stats(n, 4)] = corr(xd, yd, 'type', 'Pearson');
end

% axis
axis([xrange yrange])
set(gca, 'box', 'off', 'TickDir', 'out')
set(gca, 'YTick', yrange)
if xrange(1) < 0 && xrange(2) > 0
    hold on;
    plot([0 0], yrange, '-', 'color', 0.5*[1 1 1])
    set(gca, 'XTick', [xrange(1) 0 xrange(2)])
else
    set(gca, 'XTick', xrange)
end
if yrange(1) < 0 && yrange(2) > 0
    hold on;
    plot(xrange, [0 0], '-', 'color', 0.5*[1 1 1])
    set(gca, 'YTick', [yrange(1) 0 yrange(2)])
else
    set(gca, 'YTick', yrange)
end
% axis square

% display stats (parametric)
for n = 1:n_group
    g_lab = ['group ' num2str(n)];
    if stats(n, 4) < 0.05
        text(xrange(1)+0.025*(xrange(2)-xrange(1)), yrange(1)+(0.975 - 0.075*(n-1))*(yrange(2)-yrange(1)), ...
            [g_lab ': n=' num2str(stats(n, 5)) ', r=' num2str(stats(n, 3)) ', p<' num2str(pval_inequality(stats(n, 4)))], ...
            'fontsize', 6, 'color', map(n,:))
    else
        text(xrange(1)+0.025*(xrange(2)-xrange(1)), yrange(1)+(0.975 - 0.075*(n-1))*(yrange(2)-yrange(1)), ...
            [g_lab ': n=' num2str(stats(n, 5)) ', r=' num2str(stats(n, 3)) ', p=' num2str(pval_inequality(stats(n, 4)))], ...
            'fontsize', 6, 'color', map(n,:))
    end
end

% permutation test
if permutation_option==1 && n_group==2
    repeats = 1000;
    r_delta = nan(1, repeats);    
    x1 = x(label==groups(1)); y1 = y(label==groups(1));
    x2 = x(label==groups(2)); y2 = y(label==groups(2));
    nans = isnan(x1) | isnan(y1) | isnan(x2) | isnan(y2);
    x1(nans) = []; y1(nans) = []; x2(nans) = []; y2(nans) = [];
    lenv = length(x1);
    for r = 1:repeats
        % resample
        idxv = randi(lenv, lenv, 1);
        
        % new corr
        r1 = corr(x1(idxv), y1(idxv), 'type', 'Spearman');
        r2 = corr(x2(idxv), y2(idxv), 'type', 'Spearman');
        
        % difference
        r_delta(r) = r1 - r2;
    end        
    % compute p-value
    p_perm = [sum(r_delta < 0), sum(r_delta > 0)]./repeats;
    
    % display
    if p_perm(1) < p_perm(2)
        text(xrange(1)+0.025*(xrange(2)-xrange(1)), yrange(1)+(0.975 - 0.075*n)*(yrange(2)-yrange(1)), ...
            ['perm test: (g1 > g2) p=' num2str(p_perm(1))], 'fontsize', 6, 'color', [0 0 0])
    elseif p_perm(1) > p_perm(2)
        text(xrange(1)+0.025*(xrange(2)-xrange(1)), yrange(1)+(0.975 - 0.075*n)*(yrange(2)-yrange(1)), ...
            ['perm test: (g1 < g2) p=' num2str(p_perm(2))], 'fontsize', 6, 'color', [0 0 0])
    else
        text(xrange(1)+0.025*(xrange(2)-xrange(1)), yrange(1)+(0.975 - 0.075*n)*(yrange(2)-yrange(1)), ...
            ['perm test: (g1 = g2) p=' num2str(p_perm(2))], 'fontsize', 6, 'color', [0 0 0])
    end
end