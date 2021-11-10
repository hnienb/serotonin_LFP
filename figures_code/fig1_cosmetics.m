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
p = uipanel('Units', 'centimeters', 'Position', [xbegin+figspace_x*1.5 ybegin+figspace_y*1.5 sq*4.5 sq*4.5], 'BackgroundColor', 'white');
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
%gain_plot.XLabel.FontName = 'Arial';
%gain_plot.XLabel.FontSize = 6;
%gain_plot.YLabel.FontName = 'Arial';
%gain_plot.YLabel.FontSize = 6;

for i=1:3
   gain_plot(i).FontSize = 6;
   gain_plot(i).LineWidth = 0.5;
end

% Central scatterplot
gain_plot(3).Children(1).FontSize = 6;
gain_plot(3).Children(4).FontSize = 6;

% Right-side histogram
gain_plot(2).Children(1).FontSize = 6;
gain_plot(2).Children(3).FontSize = 6;
gain_plot(2).Children(3).Position(1) = gain_plot(2).Children(3).Position(1) - 0.015;
gain_plot(2).Children(1).Position(1) = gain_plot(2).Children(1).Position(1) + 0.015;
gain_plot(2).Children(3).Position(2) = gain_plot(2).Children(3).Position(2) + 0.7;
gain_plot(2).Children(1).Position(2) = gain_plot(2).Children(1).Position(2) + 0.7;
% This string says -0.00 (due to being a very small just below 0), so fix
% this to say 0.00
gain_plot(2).Children(1).String = '0.00';
gain_plot(2).Children(2).MarkerSize = 3;
gain_plot(2).Children(4).MarkerSize = 3;

% Top-side histogram
gain_plot(1).Children(1).FontSize = 6;
gain_plot(1).Children(3).FontSize = 6;
gain_plot(1).Children(3).Position(1) = gain_plot(1).Children(3).Position(1) - 0.04;
gain_plot(1).Children(1).Position(1) = gain_plot(1).Children(1).Position(1) + 0.04;
gain_plot(1).Children(2).MarkerSize = 3;
gain_plot(1).Children(4).MarkerSize = 3;


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
a1 = axes('parent', p, 'position', expos_1, 'fontsize', 6);
set(a1, 'ytick', get(h,'ytick'), 'yticklabel', get(h,'yticklabel'))
set(a1, 'ylim', get(h,'ylim'))
set(a1, 'ylabel', get(h,'ylabel'))
set(a1, 'ticklength', [0.02 0.02])

% Top hist - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.05;
expos_2(4) = 1/10^10;
a2 = axes('parent', p, 'position', expos_2, 'fontsize', 6);
set(a2, 'xtick', get(h,'xtick'), 'xticklabel', get(h,'xticklabel'))
set(a2, 'xlim', get(h,'xlim'))
set(a2, 'xlabel', get(h,'xlabel'))
set(a2, 'ticklength', [0.02 0.02])

% Right histogram
h = gain_plot(2);
box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Right histogram - y axis ?
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.05;
expos_1(3) = 1/10^10;
a1 = axes('parent', p, 'position', expos_1, 'fontsize', 6);
set(a1, 'ytick', get(h,'xtick'), 'yticklabel', get(h,'xticklabel'))
set(a1, 'ylim', get(h,'xlim'))
set(a1, 'ticklength', [0.02 0.02])


% Right histogram - x axis ?
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.02;
expos_2(4) = 1/10^10;
a2 = axes('parent', p, 'position', expos_2, 'fontsize', 6);
set(a2, 'xtick', [0, get(h,'ytick')], 'xticklabel', [0, get(h,'yticklabel')])
set(a2, 'xlim', get(h,'ylim'))
set(a2, 'ticklength', [0.02 0.02])


% Center scatterplot
h = gain_plot(3);
box(h, 'off');
axis(h, 'off');
base_pos = get(h, 'position');

% Center scatter - y axis
expos_1 = base_pos;
expos_1(1) = expos_1(1) - expos_1(3) * 0.02;
expos_1(3) = 1/10^10;
a1 = axes('parent', p, 'position', expos_1, 'fontsize', 6);
set(a1, 'ytick', get(h,'ytick'), 'yticklabel', get(h,'yticklabel'))
set(a1, 'ylim', get(h,'ylim'))
set(a1, 'ylabel', get(h,'ylabel'))
set(a1, 'ticklength', [0.01 0.01])

% Center scatter - x axis
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4) * 0.02;
expos_2(4) = 1/10^10;
a2 = axes('parent', p, 'position', expos_2, 'fontsize', 6, 'XScale', 'log');
set(a2, 'xtick', get(h,'xtick'), 'xticklabel', get(h,'xticklabel'))
a2.XTick([2,4]) = []; a2.XTickLabel([2,4]) = [];
set(a2, 'xlim', get(h,'xlim'))
set(a2, 'xlabel', get(h,'xlabel'))
set(a2, 'ticklength', [0.01 0.01])
%a2.Position(1) = a2.Position(1) + 0.005; % This is a correction for a weird bug that shifts the axis ever so slightly
% The axis tick at 4 should be aligned with unity line x = 4


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