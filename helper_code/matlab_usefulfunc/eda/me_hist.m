function me_hist(x, params)
%%
% generate a histogram with lines of mean and median
% INPUT:
% x ... vector
% params ... cell array containing figure parameters
%
% example: me_hist(randn(100,1),{'markerfacecolor','r','markerfacealpha',0.4})
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 2
    params = {'FaceColor','b','FaceAlpha',0.4,'EdgeColor','w',...
        'EdgeAlpha',0.8, 'linewidth',0.5};
end

% histogram
histogram(x, params{:})

% stats
xx = get(gca, 'XLim');
yy = get(gca, 'YLim');
hold on;
plot(mean(x)*[1 1], yy, ':', 'color', 0.5*[1 1 1], 'linewidth',1.5)
hold on;
plot(median(x)*[1 1], yy, '--', 'color', 0.5*[1 1 1], 'linewidth',1.5)
[~,p] = ttest(x);
text(xx(end) - 0.2*(xx(end)-xx(1)), yy(2) - 0.15*(yy(end)-yy(1)), ['p_{ttest} = ' num2str(p)])
p = signrank(x);
text(xx(end) - 0.2*(xx(end)-xx(1)), yy(2) - 0.05*(yy(end)-yy(1)), ['p_{signrank} = ' num2str(p)])

% axis
ylim(yy)
set(gca, 'YTick', yy)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')