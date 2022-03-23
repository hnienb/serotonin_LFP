function waterfallchart4mat(data, varargin)
%%
% extension of the waterfallchart for matrix data
%
%WATERFALLCHART Waterfall chart.
%   WATERFALLCHART(X) draws a waterfall chart for the values in the row
%   vector X. The chart shows the relative changes at each step between the
%   first and last elements. The colors are set by the colormap by default.
%
%   WATERFALLCHART(X,'width',W) sets the width of the bars. Value of W must
%   be between 0 and 1. The default value is W = 0.6.
%
%   h = WATERFALLCHART(...) returns a structure containing the handles to
%   the patch object (in handles.patch) and line objects (in
%   handles.lines). 
%
%   Examples: waterfallchart([1 3 4 2 5 3 4]);
%             waterfallchart(rand(1,5), 'width', 0.3);
%
%   Author: Patrick Kalita
%

% waterfall chart
me = mean(data, 1);
sem = std(data, 1)/sqrt(size(data, 1));
h = waterfallchart(me);
hold on;
errorbar(1:length(me), me, sem, 'LineStyle','none', 'color', 'k', 'capsize', 0)

