%%
% figure template 

close all

% set figure and axis
[figPars, axPars] = setPlotPars;
% figPos = [10 10 21 29.7]; % this is in cm
figPos = [5 5 21 20]; % this is in cm
figure(figPars);
set(gcf, 'position', figPos, 'paperposition', figPos);

xbegin = 2.5;
ybegin = 16;
sq = 2.5;
offset_figlab = 1.8;
figspace_x = 5;
figspace_y = 1.5;

%%
% panels 2x3 =========================
% replace 'your figure path' with your figure path
%%
% A
ax_new = axes(axPars, 'position', [xbegin ybegin sq sq]);
% open figures
fig = openfig('your figure path','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
offset_axis(0.05, axPars)

axes(axPars,'position',[xbegin-offset_figlab ybegin+sq-1 1 1])
title('A','fontsize',8)
axis off


% B
ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin sq sq]);
% open figures
fig = openfig('your figure path','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
offset_axis(0.05, axPars)

axes(axPars,'position',[xbegin-offset_figlab+figspace_x ybegin+sq-1 1 1])
title('B','fontsize',15)
axis off

% C
ax_new = axes(axPars, 'position', [xbegin+2*figspace_x ybegin sq sq]);
% open figures
fig = openfig('your figure path','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
offset_axis(0.05, axPars);

axes(axPars,'position',[xbegin+2*figspace_x-offset_figlab ybegin+sq-1 1 1])
title('C','fontsize',15)
axis off

% D
ax_new = axes(axPars, 'position', [xbegin ybegin-figspace_y sq sq]);
% open figures
fig = openfig('your figure path','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
offset_axis(0.05, axPars)

axes(axPars,'position',[xbegin-offset_figlab ybegin+sq-1-figspace_y 1 1])
title('D','fontsize',8)
axis off

% E
ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin-figspace_y sq sq]);
% open figures
fig = openfig('your figure path','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
offset_axis(0.05, axPars)

axes(axPars,'position',[xbegin-offset_figlab+figspace_x ybegin+sq-1-figspace_y 1 1])
title('E','fontsize',8)
axis off

% F
ax_new = axes(axPars, 'position', [xbegin+2*figspace_x ybegin-figspace_y sq sq]);
% open figures
fig = openfig('your figure path','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
offset_axis(0.05, axPars);

axes(axPars,'position',[xbegin+2*figspace_x-offset_figlab ybegin+sq-1-figspace_y 1 1])
title('F','fontsize',8)
axis off

%%
% save
figname = 'your figure path';
savefig([figname '.fig'])
print(gcf,'-dpdf',[figname '.pdf'],sprintf('-r%d',300))