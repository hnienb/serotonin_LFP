% Converts the Figure 1 subplots into the paper figure.

close all

% Get the 'relative' folder path to get the resource folderpath
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-2), '/');

addpath(genpath([dir_path, '/external_libraries/']));
addpath(genpath([dir_path, '/helper_code/']));


[figPars, axPars] = setPlotPars;
figPos = [5 5 21 20]; % this is in cm
figure(figPars);
set(gcf, 'position', figPos, 'paperposition', figPos);

xbegin = 2;
ybegin = 16;
sq = 2.2;
offset_figlab = 1.8;
figspace_x = 3.8;
figspace_y = -3.8;

figpath = [dir_path, '/resources/Figure_1/'];
addpath(figpath)
cols = [0 0 0; 1 0 0];
mrcs = {'o', 'square'};
msz = {7, 10};


%%
% first row, first column
%
% ==========================================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin-figspace_y/2 2*sq sq/3]);

% open figures
fig = openfig([figpath '/subplots/ka_258_20_base.fig'], 'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(1);

% copy the old object to the new axis
copyobj(ax_old.Children, ax_new); delete(fig);

% Plot a depiction of 100ms for time reference
hold on
plot([1.800, 1.900], [-200, -200], 'LineWidth', 1.25, 'Color', 'Black')
hold off

% Make axis invisible
axis off 

% ==========================================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin-figspace_y/4 2*sq sq/3]);

% open figures
fig = openfig([figpath '/subplots/ka_258_20_5HT.fig'], 'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(1);

% copy the old object to the new axis
copyobj(ax_old.Children, ax_new); delete(fig);
axis off

%%
% second row, first column
%
% ==========================================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin+figspace_y sq sq]);

% open figures
fig = openfig([figpath '/subplots/fr.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 

% extract data
data{1} = {[axesObjs(2).Children(4).XData; axesObjs(2).Children(4).YData], ...
    [axesObjs(1).Children(4).XData; axesObjs(1).Children(4).YData]};
data{2} = {[axesObjs(2).Children(5).XData; axesObjs(2).Children(5).YData], ...
    [axesObjs(1).Children(5).XData; axesObjs(1).Children(5).YData]}; 

% delete fig
delete(fig);

% plot
for a = 1:2 % animal
    for d = 1:2 % drug
        s = scatter(ax_new, data{a}{d}(1, :), data{a}{d}(2, :), msz{a}, 'marker', mrcs{a}, ...
            'markerfacecolor', cols(d, :), 'markeredgecolor', cols(d, :), 'markerfacealpha', 0.4, ...
            'markeredgealpha', 0.4, 'linewidth', 0.05);
        hold on;
    end
end
range = [0 120];
hold on;
plot(range, range, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25)

% format
xlim(range)
ylim(range)
set(gca, 'XTick', range)
set(gca, 'YTick', range)
xlabel('baseline', 'fontsize', 6)
ylabel('drug', 'fontsize', 6)
text(range(1)+0.01*(range(2)-range(1)), range(2) - 0.05*(range(2)-range(1)), ...
    ['n = ' num2str(size(data{1}{2}, 2) + size(data{2}{2}, 2))], 'color', cols(2, :), 'fontsize', 6)
text(range(1)+0.01*(range(2)-range(1)), range(2) - 0.15*(range(2)-range(1)), ...
    ['n = ' num2str(size(data{1}{1}, 2) + size(data{2}{1}, 2))], 'color', cols(1, :), 'fontsize', 6)
title({'mean firing rate', '(spk/s)'}, 'fontsize', 6)

current_ax = gca; % Obtain figure handle

% Change all text font to meet specifications
current_ax.Children(1).FontName = "Arial";
current_ax.Children(2).FontName = "Arial";
current_ax.Title.FontName = "Arial";
current_ax.FontName = "Arial";
current_ax.XLabel.FontName = "Arial";
current_ax.YLabel.FontName = "Arial";

% offset axis
offset_axis(0.05, axPars)

%%
% second row, second column
main_fig = gcf;

% open figures
fig = openfig([figpath '/subplots/gain_change.fig'], 'invisible');

% generate 'axis' (tiledlayout) object
axesObjs = get(fig, 'Children');

% Copy over each individual axis object to the main figure
center_ax = axesObjs(1);
side_ax = axesObjs(2);
top_ax = axesObjs(3);
center_ax.Parent = main_fig;
side_ax.Parent = main_fig;
top_ax.Parent = main_fig;
center_ax.Units = 'centimeters';
side_ax.Units = 'centimeters';
top_ax.Units = 'centimeters';
center_ax.Position = [xbegin+figspace_x ybegin+figspace_y sq sq];
top_ax.Position = [xbegin+figspace_x ybegin+figspace_y+sq+0.2 sq sq/3];
side_ax.Position = [xbegin+figspace_x+sq+0.2 ybegin+figspace_y sq/3 sq];

% Manually fix the visual parameters for readability
center_ax.FontSize = 6;
center_ax.FontName = "Arial";
center_ax.LineWidth = 0.5;
top_ax.FontSize = 6;
top_ax.FontName = "Arial";
top_ax.LineWidth = 0.5;
side_ax.FontSize = 6;
side_ax.FontName = "Arial";
side_ax.LineWidth = 0.5;


% Central scatterplot

% Narrow down y-range to exclude outlier point (focus on other data points)
center_ax.YLim = [-0.68, 0.68];
center_ax.YTick = [-0.68, 0, 0.68];

% Move the group number labels down on the plot after fixing font size
center_ax.Children(1).FontSize = 6;
center_ax.Children(4).FontSize = 6;
center_ax.Children(1).FontName = "Arial";
center_ax.Children(4).FontName = "Arial";
center_ax.Children(1).Position(2) = 0.60;
center_ax.Children(4).Position(2) = 0.47;
center_ax.XLabel.FontSize = 6;
center_ax.YLabel.FontSize = 6;
center_ax.XLabel.FontName = "Arial";
center_ax.YLabel.FontName = "Arial";
center_ax.XColor = [0 0 0];
center_ax.YColor = [0 0 0];
center_ax.YLabel.Color = [0 0 0];
center_ax.XLabel.Color = [0 0 0];


% Change color and thickness of unity lines to be more uniform with other figures
for i = 7:8
    center_ax.Children(i).LineWidth = 0.25;
	center_ax.Children(i).Color = [0.5, 0.5, 0.5];
end

% Reduce the size of the markers and edges to make uniform with other figures
for i = [2,3,5,6]
    center_ax.Children(i).LineWidth = 0.05;
    center_ax.Children(i).MarkerEdgeColor = center_ax.Children(i).MarkerFaceColor;
    center_ax.Children(i).MarkerEdgeAlpha = center_ax.Children(i).MarkerFaceAlpha;
end
center_ax.Children(2).SizeData = 10;
center_ax.Children(3).SizeData = 7;
center_ax.Children(5).SizeData = 10;
center_ax.Children(6).SizeData = 7;



% Right-side histogram

% Change range to follow central scatterplots y-axis narrowing
side_ax.XLim = [-0.68, 0.68];

% Update the histogram bounds and bins to follow narrow 'y'-range
for i = [5,6]
    side_ax.Children(i).BinLimits = [-0.68, 0.68];
    side_ax.Children(i).NumBins = 35;
end

% Scale the 'Counts' axis to reflect the previous changes
side_ax.YLim = [0 14];
side_ax.YTick = 14;

% Change the median labels to make uniform with other figures
for i = [1,3]
    side_ax.Children(i).FontSize = 6;
    side_ax.Children(i).FontName = "Arial";
    side_ax.Children(i).FontWeight = 'normal';
    side_ax.Children(i).Position(1) = side_ax.Children(i).Position(1) - 0.055*(i-2);
    side_ax.Children(i).Position(2) = 15;
end

% This string says -0.00 (just below 0), so fix this to say 0.00
side_ax.Children(1).String = '0.00';

% Reduce the size and position of markers for medians to fit better
for i = [2,4]
    side_ax.Children(i).MarkerSize = 1;
    side_ax.Children(i).YData = side_ax.Children(i-1).Position(2) - 2;
end


% Top-side histogram

% Change the median labels to make uniform with other figures
for i = [1,3]  
    top_ax.Children(i).FontSize = 6;
    top_ax.Children(i).FontName = "Arial";
    top_ax.Children(i).FontWeight = 'normal';
    top_ax.Children(i).Position(1) = top_ax.Children(i).Position(1) + 0.15*(2-i);
    top_ax.Children(i).Position(2) = top_ax.Children(i).Position(2) - 0.1;
end

% Reduce the marker size for medians to make uniform with other figures
top_ax.Children(2).MarkerSize = 1;
top_ax.Children(4).MarkerSize = 1;

% Manually adjust the axes since offsetAxis requires gca object

% Top histogram
h = top_ax;

box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Top hist - y axis
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.05;
expos_1(3) = 1/10^10;

a1_top = axes('parent', main_fig, 'units', 'centimeters', 'position', expos_1, 'fontsize', 6, 'TickDir', 'out');
set(a1_top, 'ytick', get(h,'ytick'), 'yticklabel', get(h,'yticklabel'), ...
    'ylim', get(h,'ylim'), 'ylabel', get(h,'ylabel'), 'ticklength', [0.02 0.02])

% Top hist - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.05;
expos_2(4) = 1/10^10;

a2_top = axes('parent', main_fig, 'units', 'centimeters', 'position', expos_2, 'fontsize', 6, 'TickDir', 'out');
set(a2_top, 'xtick', get(h,'xtick'), 'xticklabel', get(h,'xticklabel'), 'xlim', ...
    get(h,'xlim'), 'xlabel', get(h,'xlabel'), 'ticklength', [0.02 0.02])


% Right histogram
h = side_ax;

box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Right histogram - y axis
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.05;
expos_1(3) = 1/10^10;

a1_right = axes('parent', main_fig, 'units', 'centimeters', 'position', expos_1, 'fontsize', 6, 'TickDir', 'out');
set(a1_right, 'ytick', get(h,'xtick'), 'yticklabel', get(h,'xticklabel'), 'ylim', ...
    get(h,'xlim'), 'ticklength', [0.02 0.02])


% Right histogram - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.05;
expos_2(4) = 1/10^10;

a2_right = axes('parent', main_fig, 'units', 'centimeters', 'position', expos_2, 'fontsize', 6, 'TickDir', 'out');
set(a2_right, 'xtick', [0, get(h,'ytick')], 'xticklabel', [0, get(h,'yticklabel')], ...
    'xlim', get(h,'ylim'), 'ticklength', [0.02 0.02])


% Center scatterplot
h = center_ax;

box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Center scatter - y axis
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.05;
expos_1(3) = 1/10^10;

a1_center = axes('parent', main_fig, 'units', 'centimeters', 'position', expos_1, 'fontsize', 6, 'TickDir', 'out');
set(a1_center, 'ytick', get(h,'ytick'), 'yticklabel', get(h,'yticklabel'), 'ylim', ...
    get(h,'ylim'), 'ylabel', get(h,'ylabel'), 'ticklength', [0.02 0.02])

% Center scatter - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.05;
expos_2(4) = 1/10^10;

a2_center = axes('parent', main_fig, 'units', 'centimeters', 'position', expos_2, 'fontsize', 6, 'XScale', 'log', 'TickDir', 'out');
set(a2_center, 'xtick', get(h,'xtick'), 'xticklabel', get(h,'xticklabel'), 'xlim', ...
    get(h,'xlim'), 'xlabel', get(h,'xlabel'), 'ticklength', [0.02 0.02])

% Remove unnecessary labels
a2_center.XTick([2,4]) = []; a2_center.XTickLabel([2,4]) = [];
a2_center.XMinorTick = 'off';

% Remove the opened figure
delete(fig);

%% Save the figure as .fig and .pdf versions

% autosave figure
savefig([figpath '/Figure_1.fig'])

% Manually setting the pdf size to fit the plot properly
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 10]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 10]);

print(gcf,'-dpdf', [figpath '/Figure_1.pdf'])