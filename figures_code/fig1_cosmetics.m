close all

%%
% Get the 'relative' folder path to get the resource folderpath
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-2), '/');

root_path = strjoin(parts(1:end-2), '/');
addpath(genpath([root_path, '/external_libraries/']));


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
msz = 10;


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
        s = scatter(ax_new, data{a}{d}(1, :), data{a}{d}(2, :), msz, 'marker', mrcs{a}, ...
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

% offset axis
offset_axis(0.05, axPars)

%%
% second column
% generate an panel
p = uipanel('Units', 'centimeters', 'Position', [xbegin+figspace_x ybegin+figspace_y*1.120 sq*1.7 sq*1.7], 'BackgroundColor', 'white');
p.BorderType = 'None';

% open figures
fig = openfig([figpath '/subplots/gain_change.fig'], 'invisible');

% generate 'axis' (tiledlayout) object
axesObjs = get(fig, 'Children');

% Copy the 'axis' into the new panel
for i = 1:3
    axesObjs(i).Parent = p;
end


% Manually fix the visual parameters for readability
gain_plot = p.Children;

for i=1:3
   gain_plot(i).FontSize = 6;
   gain_plot(i).LineWidth = 0.5;
end

% Central scatterplot
gain_plot(3).Children(1).FontSize = 6;
gain_plot(3).Children(4).FontSize = 6;

% Reduce the size of the markers and edges to make uniform with other figures
for i = [2,3,5,6]
    gain_plot(3).Children(i).SizeData = 10;
    gain_plot(3).Children(i).LineWidth = 0.05;
    gain_plot(3).Children(i).MarkerEdgeColor = gain_plot(3).Children(i).MarkerFaceColor;
    gain_plot(3).Children(i).MarkerEdgeAlpha = gain_plot(3).Children(i).MarkerFaceAlpha;
end

% Move the group number labels down on the plot
gain_plot(3).Children(1).Position(2) = 0.23;
gain_plot(3).Children(4).Position(2) = 0.18;

% Narrow down y-range to exclude outlier point (focus on other data points)
gain_plot(3).YLim = [-0.25, 0.25];
gain_plot(3).YTick = [-0.25, 0, 0.25];

% Change color and thickness of unity lines to be more uniform with other figures
for i = 7:8
    gain_plot(3).Children(i).LineWidth = 0.25;
	gain_plot(3).Children(i).Color = [0.5, 0.5, 0.5];
end

% Right-side histogram

% Change the median labels to make uniform with other figures
for i = [1,3]
    gain_plot(2).Children(i).FontSize = 5;
    gain_plot(2).Children(i).FontWeight = 'normal';
    gain_plot(2).Children(i).Position(1) = gain_plot(2).Children(i).Position(1);
    gain_plot(2).Children(i).Position(2) = gain_plot(2).Children(i).Position(2) + 0.9;
end

% This string says -0.00 (just below 0), so fix this to say 0.00
gain_plot(2).Children(1).String = '0.00';

% Reduce the marker size for medians to fit better
gain_plot(2).Children(2).MarkerSize = 1;
gain_plot(2).Children(4).MarkerSize = 1;

% Change range to follow central scatterplots y-axis narrowing
gain_plot(2).XLim = [-0.25, 0.25];

% Update the histogram bounds and bins to follow narrow 'y'-range
for i = [5,6]
    gain_plot(2).Children(i).BinLimits = [-0.25, 0.25];
    gain_plot(2).Children(i).NumBins = 20;
end

% Top-side histogram

% Change the median labels to make uniform with other figures
for i = [1,3]
    gain_plot(1).Children(i).FontSize = 5;
    gain_plot(1).Children(i).FontWeight = 'normal';
    gain_plot(1).Children(i).Position(1) = gain_plot(1).Children(i).Position(1) + 0.12*(2-i);
    gain_plot(1).Children(i).Position(2) = gain_plot(1).Children(i).Position(2) + 0.2;
end

% Reduce the marker size for medians to make uniform with other figures
gain_plot(1).Children(2).MarkerSize = 1;
gain_plot(1).Children(4).MarkerSize = 1;

% Manually adjust the axes since offsetAxis requires gca object

% Top histogram
h = gain_plot(1);
box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Top hist - y axis
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.02;
expos_1(3) = 1/10^10;
a1_top = axes('parent', p, 'position', expos_1, 'fontsize', 6, 'TickDir', 'out');
set(a1_top, 'ytick', get(h,'ytick'), 'yticklabel', get(h,'yticklabel'), ...
    'ylim', get(h,'ylim'), 'ylabel', get(h,'ylabel'), 'ticklength', [0.02 0.02])

% Top hist - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.05;
expos_2(4) = 1/10^10;
a2_top = axes('parent', p, 'position', expos_2, 'fontsize', 6, 'TickDir', 'out');
set(a2_top, 'xtick', get(h,'xtick'), 'xticklabel', get(h,'xticklabel'), 'xlim', ...
    get(h,'xlim'), 'xlabel', get(h,'xlabel'), 'ticklength', [0.02 0.02])


% Right histogram
h = gain_plot(2);
box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Right histogram - y axis
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.05;
expos_1(3) = 1/10^10;
a1_right = axes('parent', p, 'position', expos_1, 'fontsize', 6, 'TickDir', 'out');
set(a1_right, 'ytick', get(h,'xtick'), 'yticklabel', get(h,'xticklabel'), 'ylim', ...
    get(h,'xlim'), 'ticklength', [0.02 0.02])


% Right histogram - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.02;
expos_2(4) = 1/10^10;
a2_right = axes('parent', p, 'position', expos_2, 'fontsize', 6, 'TickDir', 'out');
set(a2_right, 'xtick', [0, get(h,'ytick')], 'xticklabel', [0, get(h,'yticklabel')], ...
    'xlim', get(h,'ylim'), 'ticklength', [0.02 0.02])


% Center scatterplot
h = gain_plot(3);
box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Center scatter - y axis
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.02;
expos_1(3) = 1/10^10;
a1_center = axes('parent', p, 'position', expos_1, 'fontsize', 6, 'TickDir', 'out');
set(a1_center, 'ytick', get(h,'ytick'), 'yticklabel', get(h,'yticklabel'), 'ylim', ...
    get(h,'ylim'), 'ylabel', get(h,'ylabel'), 'ticklength', [0.02 0.02])

% Center scatter - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.02;
expos_2(4) = 1/10^10;
a2_center = axes('parent', p, 'position', expos_2, 'fontsize', 6, 'XScale', 'log', 'TickDir', 'out');
set(a2_center, 'xtick', get(h,'xtick'), 'xticklabel', get(h,'xticklabel'), 'xlim', ...
    get(h,'xlim'), 'xlabel', get(h,'xlabel'), 'ticklength', [0.02 0.02])

% Remove unnecessary labels
a2_center.XTick([2,4]) = []; a2_center.XTickLabel([2,4]) = [];
a2_center.XMinorTick = 'off';


% Fix the lengths of axes so that they are equal
a1_center.Position(4) = a2_center.Position(3);
a1_right.Position(4) = a2_center.Position(3);
h.Position([3,4]) = [a2_center.Position(3), a2_center.Position(3)];

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

print(gcf,'-dpdf', [figpath '/Figure_1.pdf'], sprintf('-r%d',300))