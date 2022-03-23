function fill_between(x,y_bottom, y_top, maincolor,transparency,varargin)
%% python-matplotlib-like function using 'fill' for shaded error bar
%
%
% written by Katsuhisa (07.02.17)
% +++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 3
        error('x, y_bottom, y_top are required as input arguments')
elseif nargin==3
        maincolor = [0 0 0];
        transparency = [];
elseif nargin==4
        transparency = [];
end

edgecolor = maincolor + (1 - maincolor)*0.55;

h = fill([x fliplr(x)],[y_bottom fliplr(y_top)],edgecolor);
set(h,'EdgeColor','none')
if ~isnan(transparency)
        alpha(transparency)
end
