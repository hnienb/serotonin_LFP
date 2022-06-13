% Converts the Figure 5 subplots into the paper figure.

close all

% Get the 'relative' folder path to get the resource folderpath
path = mfilename( 'fullpath' );

if ispc % Windows file system
    parts = strsplit(path, '\');
else
    parts = strsplit(path, '/');
end

dir_path = strjoin(parts(1:end-2), '/');
addpath(genpath([dir_path, '/helper_code']));

[figPars, axPars] = setPlotPars;
figPos = [5 5 21 20]; % this is in cm
figure(figPars);
set(gcf, 'position', figPos, 'paperposition', figPos);

xbegin = 4;
ybegin = 16;
sq = 2.2;
offset_figlab = 1.8;
figspace_x = 4.5;
figspace_y = -4;

figpath = [dir_path, '/resources/Figure_5/'];
addpath(figpath)

cols = [0 0 0; 1 0 0];
mrcs = {'o', 'square'};
xlabs = {'baseline', 'baseline', 'high FR'};
ylabs = {'NaCl', '5HT', 'low FR'};
monkey = {'kaki', 'mango', 'all'};
nses = [13, 16; 4, 37];
msz = {7,10};
stac = 0.07;
gr = 1.618;
animals = {'monkey K', 'monkey M'};

% 3 x 6 =======================

%%
% STM schematic ==========================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/ka173_example200_base.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(1);

% copy the old object to the new axis
copyobj(ax_old.Children, ax_new); delete(fig);

% Fix the text label locations
ax_new.Children(3).Position = [-64 0 0];
ax_new.Children(2).Position = [-125 -2 0];
ax_new.Children(1).Position = [-135 -3 0];

ax_new.Children(1).FontName = "Arial";
ax_new.Children(2).FontName = "Arial";
ax_new.Children(3).FontName = "Arial";

% offset axis
offset_axis(0.05, axPars)
        
% STM prediction example ==========================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin+2*figspace_x ybegin gr*sq sq]);

% open figures
fig = openfig([figpath '/subplots/ka173_example200_drug.fig'], 'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(1);

% copy the old object to the new axis
copyobj(ax_old.Children, ax_new); delete(fig);

% Fix the text label locations
ax_new.Children(3).Position = [-64 0 0];
ax_new.Children(2).Position = [-125 -2 0];
ax_new.Children(1).Position = [-135 -3 0];

ax_new.Children(1).FontName = "Arial";
ax_new.Children(2).FontName = "Arial";
ax_new.Children(3).FontName = "Arial";

% offset axis
offset_axis(0.05, axPars)

% evaluation of STM models ===============================
pcol = [[0.2549    0.6706    0.3647]; [0.2549    0.6706    0.3647]; [0.2549    0.6706    0.3647];...
    [0.2549    0.6706    0.3647]; [0.2549    0.6706    0.3647]; [0.2549    0.6706    0.3647]];

% correlation & mutual information
ranges = {[0 0.6], [0 0.2]};
xydata = cell(6, 2);
for c = 1:6
    % generate an axis
    if c <= 3
        ax_new = axes(axPars, 'position', [xbegin+figspace_x*(c-1) ybegin+figspace_y sq sq]);
        i = 1;
    elseif c >= 4
        ax_new = axes(axPars, 'position', [xbegin+figspace_x*(c-4) ybegin+figspace_y*2 sq sq]);  
        i = 2;
    end
    
    % open figures
    fig = openfig([figpath '/subplots/corr_mi.fig'],'invisible');

    % generate axis object
    axesObjs = get(fig, 'Children'); 

    % extract data
    axo = axesObjs(end - c + 1); 
    condtitle = axo.Title.String;
    lenc = length(axo.Children);
    allstats = axo.Children(lenc-7).String(6:end);
    for a = 1:2
        % store for stats
        xydata{c, a} = [axo.Children(lenc-2-a).XData; axo.Children(lenc-2-a).YData];
        
        % plot
        s = scatter(ax_new, xydata{c, a}(1, :), xydata{c, a}(2, :), msz{a}, 'marker', mrcs{a}, ...
            'markerfacecolor', pcol(c, :), 'markeredgecolor', pcol(c, :), 'markerfacealpha', 0.4, ...
            'markeredgealpha', 0.4, 'linewidth', 0.05);
        hold(ax_new, 'on');
    end
    plot(ax_new, ranges{i}, ranges{i}, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25)

    delete(fig);
    
    % format
    xlim(ranges{i})
    ylim(ranges{i})
    set(gca, 'XTick', ranges{i}, 'XTickLabel', ranges{i})
    set(gca, 'YTick', ranges{i}, 'YTickLabel', ranges{i})

    if c < 4
        xlabel(xlabs{c}, 'fontsize', 6, "fontName", "Arial")
        ylabel(ylabs{c}, 'fontsize', 6, "fontName", "Arial")
    else
        xlabel(xlabs{c-3}, 'fontsize', 6, "fontName", "Arial")
        ylabel(ylabs{c-3}, 'fontsize', 6, "fontName", "Arial")
    end
    text(ranges{i}(1)+0.01*(ranges{i}(2)-ranges{i}(1)), ranges{i}(2) - 0.05*(ranges{i}(2)-ranges{i}(1)), ...
        allstats, 'color', 'k', 'fontsize', 6, "fontName", "Arial")
    
    % correlation
    if c==2
        title('Pearson correlation coefficient', 'fontsize', 6, 'fontName', "Arial")
    elseif c==5
        title('Information gain (bits per bin)', 'fontsize', 6, 'fontName', "Arial")
    end

    % offset axis
    offset_axis(0.05, axPars)
end
delete(fig);

% display stats
disp('correlation ++++++++++++++++++++++++++++')
for a = 1:3
    if a < 3
        statsst = pair_tests([xydata{1, a}(1,:)', xydata{1, a}(2,:)'], [xydata{2, a}(1,:)', xydata{2, a}(2,:)']);
    else
        statsst = pair_tests([[xydata{1, 1}(1,:)', xydata{1, 1}(2,:)']; [xydata{1, 2}(1,:)', xydata{1, 2}(2,:)']], ...
            [[xydata{2, 1}(1,:)', xydata{2, 1}(2,:)']; [xydata{2, 2}(1,:)', xydata{2, 2}(2,:)']]);
    end
    disp([monkey{a} ' ======================'])
    disp(['base vs NaCl: p = ' num2str(statsst.pair(1).signrank.p)])
    disp(['base vs 5HT: p = ' num2str(statsst.pair(2).signrank.p)])
    disp(['NaCl vs 5HT: p = ' num2str(statsst.ranksum.p)])
end
disp('information gain ++++++++++++++++++++++++++++')
for a = 1:3
    if a < 3
        statsst = pair_tests([xydata{4, a}(1,:)', xydata{4, a}(2,:)'], [xydata{5, a}(1,:)', xydata{5, a}(2,:)']);
    else
        statsst = pair_tests([[xydata{4, 1}(1,:)', xydata{4, 1}(2,:)']; [xydata{4, 2}(1,:)', xydata{4, 2}(2,:)']], ...
            [[xydata{5, 1}(1,:)', xydata{5, 1}(2,:)']; [xydata{5, 2}(1,:)', xydata{5, 2}(2,:)']]);
    end
    disp([monkey{a} ' ======================'])
    disp(['base vs NaCl: p = ' num2str(statsst.pair(1).signrank.p)])
    disp(['base vs 5HT: p = ' num2str(statsst.pair(2).signrank.p)])
    disp(['NaCl vs 5HT: p = ' num2str(statsst.ranksum.p)])
end


%%
% autosave figure
savefig([figpath '/Figure_4.fig'])

% Manually setting the pdf size to fit the plot properly
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 10]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 10]);

print(gcf,'-dpdf', [figpath '/Figure_4.pdf'])