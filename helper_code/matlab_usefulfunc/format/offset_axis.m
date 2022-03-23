% % EXAMPLE: 
% [figPars, axPars] = setPlotPars;
% figPos = [5 5 21 10]; % this is in cm
% figure(figPars);
% set(gcf, 'position', figPos, 'paperposition', figPos);
% 
% axes(axPars, 'position', [5 5 3 3]);
% plot(rand(1,100))
% offset_axis(0.05, axPars)

function [a1, a2] = offset_axis(shift_percentage,axPars)
%% to give 'offset axis' to your figure
% INPUT: shift_percentage ... a vector of 1 or 2 elements representing a percentage of shift. When a vector of
%             2 elements was given, its first element is used for x-axis
%             and the second is for y-axis. 
%             axPars ... axis parameter given by 'setPlotPars'
% 
% OUTPUT: a1 ... new yaxis, a2 ... new xaxis
%
% EXAMPLE is given above of this funciton.
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

h = gca;
%shift_percentage = 0.2;0.05;

box off
axis off

if length(shift_percentage)==2
        xshift_percentage = shift_percentage(1);
        yshift_percentage = shift_percentage(2);
elseif length(shift_percentage)==1
        xshift_percentage = shift_percentage;
        yshift_percentage = shift_percentage;
end
        

base_pos = get(h,'position');

expos_1 = base_pos;
expos_1(1) = expos_1(1)-expos_1(3)*xshift_percentage;
% set(a,'xlim',[-110 220]);
% set(a,'ylim',[-0 2.5]);
% set(b,'xlim',[-110 220]+bkmask_delay);




expos_1(3) = 1/10^10;
a1=axes(axPars,'position', expos_1);
expos_2 = base_pos;
expos_2(2) = expos_2(2) - expos_2(4)*yshift_percentage;
expos_2(4) = 1/10^10;
a2=axes(axPars,'Position', expos_2);



set(a2,'xtick',get(h,'xtick'),'xticklabel',get(h,'xticklabel'))
set(a2,'xlim',get(h,'xlim'))
set(a2,'xlabel',get(h,'xlabel'))
set(a2,'ticklength',[0.02 0.02])

set(a1,'ytick',get(h,'ytick'),'yticklabel',get(h,'yticklabel'))
set(a1,'ylim',get(h,'ylim'))
set(a1,'ylabel',get(h,'ylabel'))
set(a1,'ticklength',[0.02 0.02])
% set(a2)