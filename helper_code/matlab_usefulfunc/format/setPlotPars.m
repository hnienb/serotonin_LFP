function [figPars, axPars] = setPlotPars(type)
% function [figPars, axPars] = setPlotPars(type)
%  sets some default parameters for making figures and axes
%
% function [figPars, axPars] = setPlotPars
%  assumes you are writing a paper
%
% function [figPars, axPars] = setPlotPars(type), where type is 'paper',
% 'poster', 'slide'


if nargin < 1
    type = 'paper';
end

figPars.Color =         'w';
figPars.Units =         'centimeters';
figPars.PaperUnits =    'centimeters';

% since 2014b 
axPars.xcolor = 'k';
axPars.ycolor = 'k'; 

axPars.TickDir =        'out';
axPars.TickLength =     [0.1 0.1];
axPars.Box =            'off';
axPars.Units = 'centimeters';
axPars.FontName=  'Arial';

switch type
    case 'paper'
        axPars.FontSize     = 6; %formerly 8
        axPars.LineWidth    = 0.5;
    case 'poster'
        axPars.FontSize     = 24;
        axPars.LineWidth    = 2;
    case 'slide'
        axPars.FontSize =       10;
        axPars.LineWidth =       1;
    case 'otherwise'
        error('<setPlotPars> Unknown plot type');
end