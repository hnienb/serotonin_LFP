function interaction_plot(mat, unity, varargin)
% scatter plot to visualize an interaction effect
% INPUT: mat ... matrix where row represents observations and
%        column must be four-dimension vector (column 1 - 4)
%        The relationship of the columns must be summarized as following:
%               group A     group B
% effect 1        1           3
% effect 2        2           4
%
% This function creates a scatter where x-axis is effect 1 and y-axis is
% effect 2.
%
% written by Katsuhisa (05.04.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin<2
    unity = 1;
end
map = parula(5);
params0 = {50, 'filled','marker','o','markerfacecolor',map(2,:),'markerfacealpha',0.4,...
    'markeredgecolor','w','markeredgealpha',0.4};
params1 = {50, 'filled','marker','^','markerfacecolor',map(4,:),'markerfacealpha',0.4,...
    'markeredgecolor','w','markeredgealpha',0.4};

% range
xminval = min(min(mat(:, [1, 3])));
xmaxval = max(max(mat(:, [1, 3])));
yminval = min(min(mat(:, [2, 4])));
ymaxval = max(max(mat(:, [2, 4])));
xdist = xmaxval - xminval;
xrange = [xminval - 0.05*xdist, xmaxval + 0.05*xdist];
ydist = ymaxval - yminval;
yrange = [yminval - 0.05*ydist, ymaxval + 0.05*ydist];

% unity
plot(xrange, [0, 0], '-','color',0.4*ones(1,3))
hold on;
plot([0, 0], yrange, '-','color',0.4*ones(1,3))
if unity==1
    hold on;
    plot(xrange, yrange, '-','color',0.4*ones(1,3))
end

% scatter
nob = size(mat, 1);
for i = 1:nob
    hold on;
    p = plot(mat(i,[1 3]), mat(i,[2 4]), '-','color',map(3,:),'linewidth',0.25);
    p.Color(4) = 0.1;
    hold on;
    scatter(mat(i,1), mat(i,2), params0{:})
    hold on;
    scatter(mat(i,3), mat(i,4), params1{:})
end

% axis
axis([xrange yrange])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(gca, 'XTick', [xrange(1) xrange(end)])
set(gca, 'YTick', [yrange(1) yrange(end)])
% axis square

% stats
xdist = xrange(end) - xrange(1);
ydist = yrange(end) - yrange(1);
% [~,p] = ttest(x,y);
text(0.6*xdist+xrange(1),0.05*ydist+yrange(1),['n=' num2str(nob)])
p = signrank(mat(:,1),mat(:,2));
text(0.6*xdist+xrange(1),0.15*ydist+yrange(1),['p1' psign(p) num2str(pval_inequality(p))])
p = signrank(mat(:,3),mat(:,4));
text(0.6*xdist+xrange(1),0.25*ydist+yrange(1),['p2' psign(p) num2str(pval_inequality(p))])
p = signrank(mean(mat(:,[1 2]), 2), mean(mat(:,[3 4]), 2));
text(0.6*xdist+xrange(1),0.35*ydist+yrange(1),['p1vs2' psign(p) num2str(pval_inequality(p))])

function s = psign(p)
if p < 0.05
    s = '<';
else
    s = '=';
end