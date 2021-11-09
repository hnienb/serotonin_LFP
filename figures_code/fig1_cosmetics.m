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
p = uipanel('Units','centimeters', 'Position', [xbegin+figspace_x*1.5 ybegin+figspace_y sq*3 sq*3], 'BackgroundColor', 'white');
p.BorderType = 'None';

% open figures
fig = openfig([figpath '/subplots/gain_change.fig'], 'invisible');

% generate 'axis' (tiledlayout) object
axesObjs = get(fig, 'Children');

% Copy the 'axis' into the new panel
axesObjs(1).Parent = p;

% Manually fix the visual parameters for readability
tile_plot = p.Children;
tile_plot.XLabel.FontName = 'Arial';
tile_plot.XLabel.FontSize = 6;
tile_plot.YLabel.FontName = 'Arial';
tile_plot.YLabel.FontSize = 6;

for i=1:3
   tile_plot.Children(i).FontSize = 6;
   tile_plot.Children(i).LineWidth = 0.5;
end

% Central scatterplot
tile_plot.Children(1).Children(1).FontSize = 6;
tile_plot.Children(1).Children(4).FontSize = 6;

% Right-side histogram
tile_plot.Children(2).Children(1).FontSize = 6;
tile_plot.Children(2).Children(3).FontSize = 6;
tile_plot.Children(2).Children(3).Position(1) = tile_plot.Children(2).Children(3).Position(1) - 0.01;
tile_plot.Children(2).Children(1).Position(1) = tile_plot.Children(2).Children(1).Position(1) + 0.01;
% This string says -0.00 (due to being a very small just below 0), so fix
% this to say 0.00
tile_plot.Children(2).Children(1).String = '0.00';
tile_plot.Children(2).Children(2).MarkerSize = 3;
tile_plot.Children(2).Children(4).MarkerSize = 3;

% Top-side histogram
tile_plot.Children(3).Children(1).FontSize = 6;
tile_plot.Children(3).Children(3).FontSize = 6;
tile_plot.Children(3).Children(3).Position(1) = tile_plot.Children(3).Children(3).Position(1) - 0.02;
tile_plot.Children(3).Children(1).Position(1) = tile_plot.Children(3).Children(1).Position(1) + 0.02;
tile_plot.Children(3).Children(2).MarkerSize = 3;
tile_plot.Children(3).Children(4).MarkerSize = 3;

delete(fig);

%% Save the figure as .fig and .pdf versions

% autosave figure
savefig([figpath '/Figure_1.fig'])

% Manually setting the pdf size to fit the plot properly
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 7.5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 7.5]);

print(gcf,'-dpdf', [figpath '/Figure_1.pdf'], sprintf('-r%d',300))