function handles = waterfallchart(data,varargin)
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
% @ Katsuhisa modified axis range and color
%


p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && size(x,1) == 1);
p.addParamValue('width', 0.6, @(x) isnumeric(x) && x > 0 && x < 1);
p.parse(data, varargin{:});

cla;

width = p.Results.width;
data = [0 p.Results.data];
[numRows, numCols] = size(data);

k = [-width/2; -width/2; width/2; width/2];
x = meshgrid(1:numCols-1, 1:4);
x = bsxfun(@plus, x, k);

ia = 1:numCols-1;
ib = 2:numCols;
y = data([ia;ib;ib;ia]);

% rgb = 0.5*eye(3);
mycol = [0.7020    0.8039    0.8902; ...
    0.9843    0.7059    0.6824; ...
    0.8000    0.9216    0.7725];
ind = (diff(data) > 0) + 1;
ind([1 end]) = [3 3];
colormap(mycol);
handles.patch = patch(x,y,ind);
set(handles.patch, 'CDataMapping', 'direct');
xl = [x(4:4:end-1); x(5:4:end)];
yl = data([1,1],2:end-1);
handles.lines = line(xl,yl,'color','black', 'linewidth', 0.25);
end