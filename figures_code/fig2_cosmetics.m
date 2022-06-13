% Converts the Figure 2 subplots into the paper figure.

close all

% Get the 'relative' folder path to get the resource folderpath
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-2), '/');
addpath(genpath([dir_path, '/helper_code/']));

[figPars, axPars] = setPlotPars;
figPos = [5 5 21 20]; % this is in cm
figure(figPars);
set(gcf, 'position', figPos, 'paperposition', figPos);

xbegin = 1;
ybegin = 17;
sq = 2.5;
offset_figlab = 1.8;
figspace_x = 5.5;
figspace_y = -3.5;

figpath = [dir_path, '/resources/Figure_2/'];
addpath(figpath)

cols = [0 0 0; 1 0 0];
mrcs = {'o', 'square'};
drugs = {'NaCl', '5HT'};
monkey = {'monkey K', 'monkey M'};
nses = [13, 16; 4, 35];
msz = {7, 10};
map = jet;

gr = 1.618;
sq = sq/sqrt(gr);
crange = [-0.008 0.008];

%%
% colormap ==========================================================
% NaCl =======================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/spectrogram.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(4);

% get data
x = ax_old.Children(2).XData;
y = ax_old.Children(2).YData;
c = ax_old.Children(2).CData;
delete(fig);

% plot
hold(ax_new, 'on')
h = pcolor(ax_new, x, y, c);
h.EdgeColor = 'none';

% format
caxis(crange)
colormap(map)
hold on
plot(ax_new, [0 0], [3 48], '--', 'linewidth', 0.5, 'Color', [0.5 0.5 0.5]);
xlim([-0.1 2])
ylim([3 48])
set(gca, 'XTick', [0 0.8 2])
set(gca, 'YTick', [3 10 48])
set(gca, 'YScale', 'log')
set(gca, 'FontName', 'Arial')
title('baseline - NaCl (n = 17)', 'fontsize', 6, 'fontname', 'Arial')
xlabel('time after stimulus onset (sec)', 'fontsize', 6, 'fontname', "Arial")
ylabel('frequency (Hz)', 'fontsize', 6, 'fontname', "Arial")

% offset axis
a1 = offset_axis(0.05, axPars);
set(a1, 'YScale', 'log')
set(a1, 'YMinorTick', 'off')



% ==========================================================
% 5HT =======================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin+figspace_y gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/spectrogram.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(1);

% get data
x = ax_old.Children(2).XData;
y = ax_old.Children(2).YData;
c = ax_old.Children(2).CData;
delete(fig);

% plot
h = pcolor(ax_new, x, y, c);
h.EdgeColor = 'none';
hold(ax_new, 'on')

% format
caxis(crange)
c = colorbar('eastoutside');
c.Label.String = '\Delta power (a.u.)';
c.Ticks = [crange(1) 0 crange(2)];
c.TickDirection = 'out';
c.Position(1) = c.Position(1)*1.4;
c.FontSize = 6;
c.FontName = "Arial";

colormap(map)

caxis(crange)
colormap(jet)

hold on
plot(ax_new, [0 0], [3 48], '--', 'linewidth', 0.5, 'Color', [0.5 0.5 0.5]);

xlim([-0.1 2])
ylim([3 48])
set(gca, 'XTick', [0 2])
set(gca, 'YTick', [3 10 48])
set(gca, 'FontName', "Arial")
set(ax_new, 'YScale', 'log')
title('baseline - 5HT (n = 50)', 'fontsize', 6, 'fontname', 'Arial')
xlabel('time after stimulus onset (sec)', 'fontsize', 6, 'fontname', 'Arial')
ylabel('frequency (Hz)', 'fontsize', 6, 'fontname', 'Arial')

% offset axis
a1 = offset_axis(0.05, axPars);
set(a1, 'YScale', 'log')
set(a1, 'YMinorTick', 'off')


% ==========================================================
% FR =======================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin+2*figspace_y gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/spectrogram_sc.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(1);

% get data
x = ax_old.Children(2).XData;
y = ax_old.Children(2).YData;
c = ax_old.Children(2).CData;
delete(fig);

% plot
h = pcolor(ax_new, x, y, c);
h.EdgeColor = 'none';
hold(ax_new, 'on')

% format
caxis(crange)
colormap(map)

hold on
plot(ax_new, [0 0], [3 48], '--', 'linewidth', 0.5, 'Color', [0.5 0.5 0.5]);

xlim([-0.1 2])
ylim([3 48])
set(gca, 'XTick', [0 2])
set(gca, 'YTick', [3 10 48])
set(gca, 'YScale', 'log')
set(gca, 'FontName', 'Arial')
title('high FR - low FR (n = 68)', 'fontsize', 6, 'fontname', 'Arial')
xlabel('time after stimulus onset (sec)', 'fontsize', 6, 'fontname', 'Arial')
ylabel('frequency (Hz)', 'fontsize', 6, 'fontname', 'Arial')

% offset axis
a1 = offset_axis(0.05, axPars);
set(a1, 'YScale', 'log')
set(a1, 'YMinorTick', 'off')

%%
% spectrum ==========================================================
% NaCl =======================================
yy = [0 0.03];

% generate an axis
ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/psd.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(4);

% copy the old object to the new axis
copyobj(ax_old.Children, ax_new); delete(fig);

% Add line depictions for the bandranges analyzed
hold on
plot([3 10], [0 0], 'LineWidth', 2, 'Color', [0.10 0.52 1 0.4]) % Low freq
plot([30 48], [0 0], 'LineWidth', 2, 'Color', [1 0.5 0.04 0.4]) % Gamma

% format
xlim([3 48])
ylim(yy)
set(gca, 'XTick', [3 48])
set(gca, 'XScale', 'log')
set(gca, 'YTick', yy)
set(gca, 'FontName', 'Arial')
text(35, 0.024,'baseline', 'color', 'k', 'fontsize', 6, 'fontname', 'Arial')
text(35, 0.02,'NaCl', 'color', 'r', 'fontsize', 6, 'fontname', 'Arial')
xlabel('frequency (Hz)', 'fontsize', 6, 'fontname', 'Arial')
ylabel('power (a.u.)', 'fontsize', 6, 'fontname', 'Arial')

% offset axis
offset_axis(0.05, axPars);

% ==========================================================
% 5HT =======================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin+figspace_y gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/psd.fig'], 'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(2);

% copy the old object to the new axis
copyobj(ax_old.Children, ax_new); delete(fig);

% Add line depictions for the bandranges analyzed
hold on
plot([3 10], [0 0], 'LineWidth', 2, 'Color', [0.10 0.52 1 0.4]) % Low freq
plot([30 48], [0 0], 'LineWidth', 2, 'Color', [1 0.5 0.04 0.4]) % Gamma

% format
xlim([3 48])
ylim(yy)
set(gca, 'XTick', [3 48])
set(gca, 'YTick', yy)
set(gca, 'XScale', 'log')
set(gca, 'FontName', 'Arial')
text(35, 0.024,'baseline', 'color', 'k', 'fontsize', 6, 'fontname', 'Arial')
text(35, 0.02,'5HT', 'color', 'r', 'fontsize', 6, 'fontname', 'Arial')
xlabel('frequency (Hz)', 'fontsize', 6, 'fontname', 'Arial')
ylabel('power (a.u.)', 'fontsize', 6, 'fontname', 'Arial')

% offset axis
offset_axis(0.05, axPars)

% ==========================================================
% FR =======================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin+2*figspace_y gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/psd_sc.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(2);

% copy the old object to the new axis
copyobj(ax_old.Children, ax_new); delete(fig);

% Add line depictions for the bandranges analyzed
hold on
plot([3 10], [0 0], 'LineWidth', 2, 'Color', [0.10 0.52 1 0.4]) % Low freq
plot([30 48], [0 0], 'LineWidth', 2, 'Color', [1 0.5 0.04 0.4]) % Gamma

% format
xlim([3 48])
ylim(yy)

set(gca, 'XScale', 'log')
set(gca, 'XTick', [3 48])
set(gca, 'YTick', yy)
set(gca, 'FontName', "Arial")
text(35, 0.024,'high FR', 'color', 'k', 'fontsize', 6, 'fontname', 'Arial')
text(35, 0.02,'low FR', 'color', 'r', 'fontsize', 6, 'fontname', 'Arial')
xlabel('frequency (Hz)', 'fontsize', 6, 'fontname', 'Arial')
ylabel('power (a.u.)', 'fontsize', 6, 'fontname', 'Arial')

% offset axis
offset_axis(0.05, axPars)


%%
% power scatter ===============================
% ==========================================================
% drug
pcol = [0.2549    0.6706    0.3647];
fnames = {'lowfreq_scatter', 'gamma_scatter'};
tnames = {'low-frequency (\leq 10Hz)', 'gamma (\geq 30Hz)'};
tcolors = {[0.10 0.52 1], [1 0.5 0.04]};
ranges = {[0.015, 0.025; 0.005, 0.035], [0.5, 2.7; 0.1, 3.5].*0.001} ;
pvals = cell(1, 2);
for f = 1:2
    % open figures
    fig = openfig([figpath '/subplots/' fnames{f} '.fig'], 'invisible');

    % generate axis object
    axesObjs = get(fig, 'Children'); 

    % extract data
    pvals{f, 1} = axesObjs(2).Children(1).String(6:end);
    pvals{f, 2} = axesObjs(1).Children(1).String(6:end);
    data{1} = {[axesObjs(2).Children(4).XData; axesObjs(2).Children(4).YData], ...
        [axesObjs(1).Children(4).XData; axesObjs(1).Children(4).YData]};
    data{2} = {[axesObjs(2).Children(5).XData; axesObjs(2).Children(5).YData], ...
        [axesObjs(1).Children(5).XData; axesObjs(1).Children(5).YData]}; 

    % delete fig
    delete(fig);

    % plot
    for d = 1:2 % drug      
        % generate an axis
        ax_new = axes(axPars, 'position', [xbegin+(1.7+f)*figspace_x/1.4 ybegin+(d-1)*figspace_y sq sq]);
        for a = 1:2 % animal     
            % plot
            s = scatter(ax_new, data{a}{d}(1, :), data{a}{d}(2, :), msz{a}, 'marker', mrcs{a}, ...
                'markerfacecolor', pcol, 'markeredgecolor', pcol, 'markerfacealpha', 0.4, ...
                'markeredgealpha', 0.4, 'linewidth', 0.05);
            hold(ax_new, 'on');
        end
        range = 10*log10(ranges{f}(d, :));
        hold(ax_new, 'on');
        plot(ax_new, range, range, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25);

        % format
        xlim(range)
        ylim(range)
        set(gca, 'XTick', range, 'XTickLabel', ranges{f}(d, :))
        set(gca, 'YTick', range, 'YTickLabel', ranges{f}(d, :))
        set(gca, 'FontName', "Arial")
        xlabel('baseline', 'fontsize', 6, 'fontname', 'Arial')
        ylabel('drug', 'fontsize', 6, 'fontname', 'Arial')
        text(range(1)+0.01*(range(2)-range(1)), range(2) - 0.05*(range(2)-range(1)), ...
            pvals{f, d}, 'color', 'k', 'fontsize', 6, 'fontname', 'Arial')
        
        if d==1 
            title({tnames{f}, 'power (a.u.)'}, 'fontsize', 6, 'Color', tcolors{f}, 'fontname', 'Arial')
        end

        % offset axis
        offset_axis(0.05, axPars)
    end    
end

% FR
fnames = {'lowfreq_scatter_sc', 'gamma_scatter_sc'};
ranges = {[0.01, 0.03], [0.1, 5.5].*0.001} ;
pvals = cell(1,2);
for f = 1:2
    % open figures
    fig = openfig([figpath '/subplots/' fnames{f} '.fig'], 'invisible');

    % generate axis object
    axesObjs = get(fig, 'Children'); 

    % extract data
    pvals{f} = axesObjs(1).Children(1).String(6:end);
    data{1} = {[axesObjs(1).Children(4).XData; axesObjs(1).Children(4).YData]};
    data{2} = {[axesObjs(1).Children(5).XData; axesObjs(1).Children(5).YData]}; 

    % delete fig
    delete(fig);

    % plot
    % generate an axis
    ax_new = axes(axPars, 'position', [xbegin+(1.7+f)*figspace_x/1.4 ybegin+2*figspace_y sq sq]);
    for a = 1:2 % animal     
        % plot
        s = scatter(ax_new, data{a}{1}(1, :), data{a}{1}(2, :), msz{a}, 'marker', mrcs{a}, ...
            'markerfacecolor', pcol, 'markeredgecolor', pcol, 'markerfacealpha', 0.4, ...
            'markeredgealpha', 0.4, 'linewidth', 0.05);
        hold(ax_new, 'on');
    end
    range = 10*log10(ranges{f});
    hold(ax_new, 'on');
    plot(ax_new, range, range, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25);

    % format
    xlim(range)
    ylim(range)
    set(gca, 'XTick', range, 'XTickLabel', ranges{f})
    set(gca, 'YTick', range, 'YTickLabel', ranges{f})
    set(gca, 'FontName', 'Arial')
    xlabel('baseline', 'fontsize', 6, 'fontname', 'Arial')
    ylabel('drug', 'fontsize', 6, 'fontname', 'Arial')
    text(range(1)+0.01*(range(2)-range(1)), range(2) - 0.05*(range(2)-range(1)), ...
        pvals{f}, 'color', 'k', 'fontsize', 6, 'fontname', 'Arial')

    % offset axis
    offset_axis(0.05, axPars)
end    


% autosave figure
savefig([figpath '/Figure_2.fig'])

% Manually setting the pdf size to fit the plot properly
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 10]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 10]);

print(gcf,'-dpdf', [figpath '/Figure_2.pdf'])