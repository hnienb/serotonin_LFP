% Converts the Figure 3 subplots into the paper figure.

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

xbegin = 2;
ybegin = 16;
sq = 2.2;
offset_figlab = 1.8;
figspace_x = 3.8;
figspace_y = -3.8;

map = 'jet';

figpath = [dir_path, '/resources/Figure_3/'];
addpath(figpath)

cols = [0 0 0; 1 0 0];
mrcs = {'o', 'square'};
drugs = {'NaCl', '5HT'};
monkey = {'monkey K', 'monkey M'};
msz = {7,10};
nses = [13, 16; 4, 35];
stac = 0.07;
gr = 1.618;

% 3 x 6 =======================

%%
% first row
%
% STA example=================================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin ybegin sq sq]);

% open figures
fig = openfig([figpath '/subplots/sta_batch.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 
ax_old = axesObjs(1);
nspk = ax_old.Children(2).String;

% copy the old object to the new axis
copyobj(ax_old.Children(5), ax_new); delete(fig);
 
xlim(stac*[-1 1])
ylim([-0.80 0.80]) % Widened the range to properly fit the plot
set(gca, 'XTick', [-stac 0 stac])
set(gca, 'YTick', [-0.80 0 0.80])
set(gca, 'FontName', 'Arial')
xlabel('time after spike (sec)', 'fontsize', 6, 'fontname', 'Arial')
ylabel('LFP (a.u.)', 'fontsize', 6, 'fontname', 'Arial')
text(-0.065, 0.65, 'k258', 'fontsize', 6, 'fontname', 'Arial')
text(-0.065, 0.5, [num2str(nspk) ' spikes'], 'fontsize', 6, 'fontname', 'Arial')
% offset axis
offset_axis(0.05, axPars)

% stLFP scatter ==============================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin sq sq]);

% open figures
fig = openfig([figpath '/subplots/sta_amp.fig'], 'invisible');

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
range = [0 2];
hold on;
plot(range, range, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25)

% format
xlim(range)
ylim(range)
set(gca, 'XTick', range)
set(gca, 'YTick', range)
set(gca, 'FontName', "Arial")
xlabel('baseline', 'fontsize', 6, 'fontname', 'Arial')
ylabel('drug', 'fontsize', 6, 'fontname', 'Arial')
text(range(1)+0.01*(range(2)-range(1)), range(2) - 0.05*(range(2)-range(1)), ...
    ['n = ' num2str(size(data{1}{2}, 2) + size(data{2}{2}, 2))], 'color', ...
    cols(2, :), 'fontsize', 6, 'fontname', 'Arial')
text(range(1)+0.01*(range(2)-range(1)), range(2) - 0.15*(range(2)-range(1)), ...
    ['n = ' num2str(size(data{1}{1}, 2) + size(data{2}{1}, 2))], 'color', ...
    cols(1, :), 'fontsize', 6, 'fontname', 'Arial')
title('stLFP amplitude (a.u.)', 'fontsize', 6, 'fontname', 'Arial')

% offset axis
offset_axis(0.05, axPars)

% coherence ===================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin+2*figspace_x ybegin sq sq]);

% open figures
fig = openfig([figpath '/subplots/coh_scatter.fig'],'invisible');

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
xrange = [0 0.5];
yrange = [0 0.5];
hold on;
plot(xrange, yrange, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25)

% format
xlim(xrange)
ylim(yrange)
set(gca, 'XTick', [xrange(1) xrange(2)])
set(gca, 'YTick', [yrange(1) yrange(2)])
set(gca, 'FontName', 'Arial')
xlabel('baseline', 'fontsize', 6, 'fontname', 'Arial')
ylabel('drug', 'fontsize', 6, 'fontname', 'Arial')
text(xrange(1)+0.01*(xrange(2)-xrange(1)), yrange(2) - 0.05*(yrange(2)-yrange(1)), ...
    ['n = ' num2str(size(data{1}{2}, 2) + size(data{2}{2}, 2))], 'color', ...
    cols(2, :), 'fontsize', 6, 'fontname', 'Arial')
text(xrange(1)+0.01*(xrange(2)-xrange(1)), yrange(2) - 0.15*(yrange(2)-yrange(1)), ...
    ['n = ' num2str(size(data{1}{1}, 2) + size(data{2}{1}, 2))], 'color', ...
    cols(1, :), 'fontsize', 6, 'fontname', 'Arial')
title({'spike-LFP coherence', '(\leq 10Hz)'}, 'fontsize', 6, 'fontname', 'Arial')

% offset axis
offset_axis(0.05, axPars)

% delete fig
delete(fig);


% stLFP vs LFP power ===================================
% generate an axis
ax_new = axes(axPars, 'position', [xbegin+3*figspace_x ybegin sq sq]);

% open figures
fig = openfig([figpath '/subplots/sta_lowfreqpow.fig'],'invisible');

% generate axis object
axesObjs = get(fig, 'Children'); 

% extract data
data{1} = {[axesObjs(2).Children(6).XData; axesObjs(2).Children(6).YData], ...
    [axesObjs(1).Children(6).XData; axesObjs(1).Children(6).YData]};
data{2} = {[axesObjs(2).Children(8).XData; axesObjs(2).Children(8).YData], ...
    [axesObjs(1).Children(8).XData; axesObjs(1).Children(8).YData]}; 

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
xrange = [-0.5 0.6];
yrange = [-1.75 4.75]; % Fixed the yrange to fit the plot properly
hold on;
plot([0 0], yrange, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25)
hold on;
plot(xrange, [0 0], '-', 'color', 0.5*[1 1 1], 'linewidth', 0.25)

% regression line
for d = 1:2
   X = [data{1}{d}(1, :)'; data{2}{d}(1, :)'];
   beta = glmfit(X, [data{1}{d}(2, :)'; data{2}{d}(2, :)']); 
   hold on;
   plot(xrange, glmval(beta, xrange, 'identity'), '-', 'color', cols(d, :), 'linewidth', 0.5)
end

% format
xlim(xrange)
ylim(yrange)
set(gca, 'XTick', [xrange(1) 0 xrange(2)])
set(gca, 'YTick', [yrange(1) 0 yrange(2)])
set(gca, 'FontName', 'Arial')
xlabel('\Delta stLFP amplitude', 'fontsize', 6, 'fontname', 'Arial')
ylabel('\Delta low-freq power', 'fontsize', 6, 'fontname', 'Arial')
text(xrange(1)+0.01*(xrange(2)-xrange(1)), yrange(2) - 0.05*(yrange(2)-yrange(1)), ...
    ['n = ' num2str(size(data{1}{2}, 2) + size(data{2}{2}, 2))], 'color', ...
    cols(2, :), 'fontsize', 6, 'fontname', 'Arial')
text(xrange(1)+0.01*(xrange(2)-xrange(1)), yrange(2) - 0.15*(yrange(2)-yrange(1)), ...
    ['n = ' num2str(size(data{1}{1}, 2) + size(data{2}{1}, 2))], 'color', ...
    cols(1, :), 'fontsize', 6, 'fontname', 'Arial')

% offset axis
offset_axis(0.05, axPars)


%%
% stLFP and power
%
% ==========================================================
yrange = [-0.3, 0.3; -0.6 0.4]; % Changed the yranges to fit the plots
crange = [-0.001, 0.001;-0.001, 0.001];
for a = 1:2 % animals
    for d = 1:2 % drugs
        % STA =======================================
        % generate an axis
        ax_new = axes(axPars, 'position', [xbegin+figspace_x*(2*(a-1))*1.05 ybegin+figspace_y*d sq sq]);

        % open figures
        fig = openfig([figpath '/subplots/sta_rc.fig'],'invisible');

        % generate axis object
        axesObjs = get(fig, 'Children'); 
        ax_old = axesObjs(30-5*(d-1)-10*(a-1));

        % no errorbar
        for i = 1:length(ax_old.Children)
           try
               ax_old.Children(i).FaceAlpha = 0.1;
           catch
               continue
           end
        end
        
        % copy the old object to the new axis
        copyobj(ax_old.Children, ax_new); delete(fig);

        xlim(stac*[-1 1])
        ylim(yrange(a,:))
  
        set(gca, 'XTick', [-stac 0 stac])
        set(gca, 'YTick', [yrange(a,1) 0 yrange(a,2)])
        set(gca, 'FontName', 'Arial')
        xlabel('time after spike (sec)', 'fontsize', 6, 'fontname', 'Arial')
        ylabel({'LFP (a.u.)'}, 'fontsize', 6, 'fontname', 'Arial')
        
        % Fix the dashed line length
        ax_temp = gca;
        ax_temp.Children(1).YData = yrange(a,:);
        
        text(0.025, yrange(a,1)+0.95*(yrange(a,2)-yrange(a,1)), ...
            ['n = ' num2str(nses(a, d))], 'fontsize', 6, 'fontname', 'Arial')
        text(0.025, yrange(a,1)+0.85*(yrange(a,2)-yrange(a,1)), 'baseline', ...
        'color', 'k', 'fontsize', 6, 'fontname', 'Arial')
        text(0.025, yrange(a,1)+0.75*(yrange(a,2)-yrange(a,1)), drugs{d}, ...
            'color', 'r', 'fontsize', 6, 'fontname', 'Arial')
        if d==1
            title(monkey{a}, 'fontsize', 6, 'fontname', 'Arial')
        end

        % offset axis
        offset_axis(0.05, axPars)
        
        % d power =======================================
        % generate an axis
        ax_new = axes(axPars, 'position', [xbegin+figspace_x*(2*a-1) ybegin+figspace_y*d sq sq]);

        % open figures
        fig = openfig([figpath '/subplots/sta_rc.fig'],'invisible');

        % generate axis object
        axesObjs = get(fig, 'Children'); 
        ax_old = axesObjs(27-5*(d-1)-10*(a-1));

        % copy the old object to the new axis
        copyobj(ax_old.Children, ax_new); delete(fig);

        if d==1
            c = colorbar('eastoutside');
            cpos = c.Position;
            if a==1
                cpos(1) = 1.15*cpos(1);
            else
                cpos(1) = 1.08*cpos(1);
            end
            cpos(3) = 1*cpos(3);
            c.Position = cpos;
            c.AxisLocation = 'out';
            c.Box = 'off';
            c.FontSize = 6;
            c.FontName = 'Arial';
            c.TickDirection = 'out';
            c.Label.String = 'a.u.';
            c.Ticks = [crange(a,1) 0 crange(a,2)];
        end
        
        colormap(map)
        caxis(crange(a,:));
        xlim(stac*[-1 1])
        ylim([3 48])
        
        set(gca, 'XTick', [-stac 0 stac])
        set(gca, 'YTick', [3 48])
        set(gca, 'FontName', "Arial")
        xlabel('time after spike (sec)', 'fontsize', 6, 'fontname', 'Arial')
        ylabel('frequency (Hz)', 'fontsize', 6, 'fontname', 'Arial')
        title(['\Delta PSD (baseline - ' drugs{d} ')'], 'fontsize', 6, 'fontname', 'Arial')
        
        % Change color of the dashed line at 0 to black=
        child = get(gca, 'Children');
        set(child(1), 'Color', [0 0 0]);
        set(child(1), 'LineWidth', 0.5);

        % offset axis
        offset_axis(0.05, axPars)
    end
    
    % STA =======================================
    % generate an axis
    ax_new = axes(axPars, 'position', [xbegin+figspace_x*(2*(a-1))*1.05 ybegin+figspace_y*3 sq sq]);

    % open figures
    fig = openfig([figpath '/subplots/sta_sc.fig'],'invisible');

    % generate axis object
    axesObjs = get(fig, 'Children'); 
    ax_old = axesObjs(15-5*(a-1));

    % no errorbar
    for i = 1:length(ax_old.Children)
       try
           ax_old.Children(i).FaceAlpha = 0.1;
       catch
           continue
       end
    end

    % copy the old object to the new axis
    copyobj(ax_old.Children, ax_new); delete(fig);

    xlim(stac*[-1 1])
    ylim(yrange(a,:))
    set(gca, 'XTick', [-stac 0 stac])
    set(gca, 'YTick', [yrange(a,1) 0 yrange(a,2)])
    set(gca, 'FontName', "Arial")
    xlabel('time after spike (sec)', 'fontsize', 6, 'fontname', 'Arial')
    ylabel({'LFP (a.u.)'}, 'fontsize', 6, 'fontname', 'Arial')

    % Fix the dashed line length
    ax_temp = gca;
    ax_temp.Children(1).YData = yrange(a,:);
    
    % Changed the x-val on the following 3 lines to move text to the right
    text(0.025, yrange(a,1)+0.95*(yrange(a,2)-yrange(a,1)), ['n = ' num2str(sum(nses(a, :)))], 'fontsize', 6, 'fontname', 'Arial')
    text(0.025, yrange(a,1)+0.85*(yrange(a,2)-yrange(a,1)), 'high FR', 'color', 'k', 'fontsize', 6, 'fontname', 'Arial')
    text(0.025, yrange(a,1)+0.75*(yrange(a,2)-yrange(a,1)), 'low FR', 'color', 'r', 'fontsize', 6, 'fontname', 'Arial')
    
    % offset axis
    offset_axis(0.05, axPars)
    
    % d power =======================================
    % generate an axis
    ax_new = axes(axPars, 'position', [xbegin+figspace_x*(2*a-1) ybegin+figspace_y*3 sq sq]);

    % open figures
    fig = openfig([figpath '/subplots/sta_sc.fig'], 'invisible');

    % generate axis object
    axesObjs = get(fig, 'Children'); 
    ax_old = axesObjs(12-5*(a-1));

    % copy the old object to the new axis
    copyobj(ax_old.Children, ax_new); delete(fig);
    
    colormap(map)
    caxis(crange(a, :))
    xlim(stac*[-1 1])
    ylim([3 48])
 
    set(gca, 'XTick', [-stac 0 stac])
    set(gca, 'YTick', [3 48])
    set(gca, 'FontName', 'Arial')
    xlabel('time after spike (sec)', 'fontsize', 6, 'fontname', 'Arial')
    ylabel('frequency (Hz)', 'fontsize', 6, 'fontname', 'Arial')
    title(['\Delta PSD (high - low FR)'], 'fontsize', 6, 'fontname', 'Arial')
    
    % Change the dashed line at x = 0 to black
    child = get(gca, 'Children');
    set(child(1), 'Color', [0 0 0]);
    set(child(1), 'LineWidth', 0.5);

    % offset axis
    offset_axis(0.05, axPars)

end


% autosave figure
savefig([figpath '/Figure_3.fig'])

% Manually setting the pdf size to fit the plot properly
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 10]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 10]);

print(gcf,'-dpdf', [figpath '\Figure_3.pdf'])