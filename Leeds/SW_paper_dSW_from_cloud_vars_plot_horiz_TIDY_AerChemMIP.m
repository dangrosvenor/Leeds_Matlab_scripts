%% Tidy plot
yoff=0.28;
    %Changing the gca position here is the key to making the plots more
    %compact.
%set(gca,'position',[0.1367    0.1100+yoff    0.7683    0.8150-(yoff*1.2)]
%set(gca,'position',[0.3300    0.1100    0.7750/3    0.2150]);
set(gca,'position',[0.3300    0.1100    0.7750/3    0.8150]);
set(gca,'position',[0.3300    0.1100    0.7750/3    0.4150]);
h_xlab = xlabel(ylabs{ivar});
xlab_pos = get(h_xlab,'position');
%Move the xlabel down a bit to avoid is hitting the numbers on the x-axis
%set(h_xlab,'position',[xlab_pos(1) xlab_pos(2)-0.2 xlab_pos(3)]);
%.... seems to make it disappear, though...

fontsize_figure(gcf,gca,16);
title(period_lab); %set in ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot.m
grid off

x=0;
for i=2:length(section_width)
    x = x + section_width{i}; %total plot height
end
set(gca,'ylim',[0 x]);

%set(gca,'xtick',[1:length(ytick_labs)]);
set(gca,'ytick',xx);
set(gca,'yticklabel',[]);

% Add tick labels at 45 degree angles.
ylims=get(gca,'ylim');
xlims=get(gca,'xlim');

%In cases where the right x-axis limits is close to zero extend a little beyond
%zero since it looks better.
 minval = (xlims(2)-xlims(1))/20; %Min abs value of xlim
% if abs(xlims(1))<minval & sign(xlims(1))~=sign(xlims(2))
%     xlims(1) = -sign(xlims(2))*minval; %sign(0)=0 so have to deal with case when value is zero
% end
if abs(xlims(2))<minval & sign(xlims(1))~=sign(xlims(2))
    xlims(2) = -sign(xlims(1))*minval;
end
dxlim = xlims(2)-xlims(1);
xpos=xlims(1)-dxlim*0.05;
for i=1:length(ytick_labs)
    %text(xpos,xx(i),ytick_labs{i},'rotation',-45,'fontsize',14);
    text(xpos,xx(i),ytick_labs{i},'rotation',0,'fontsize',fsize_text,'HorizontalAlignment','right');
end

%Plot the horizontal line breaks
x=0;
%plot([xlims(1) xlims(2)],[x x],'k-');
for i=2:length(section_width)
    x = x + section_width{i};
    plot([xlims(1) xlims(2)],[x x],'k-');
    
    %hL=line([i-0.5 i-0.5],[ylims(1) ylims(2)]);
    %set(hL,'color','k');
end

%Plot the zero line at y=0
plot([0,0],[ylims(1) ylims(2) ],'k--');
set(gca,'YminorTick','off');
%hL=line([xlims(1) xlims(2)],[0 0]);
%set(hL,'color','k');
%

% for i=1:length(big_section_break)
%     x = big_section_break(i);
%     plot([xlims(1)-dxlim*0.5 xlims(2)],[x x],'k--','linewidth',6);
%     %hL=line([xlims(1)-dxlim*0.5 xlims(2)],[x x]);
%     %hL=line([i-0.5 i-0.5],[ylims(1) ylims(2)]);
%     %set(hL,'color','k','linestyle','--','linewidth',6);
% end

set(gca,'xlim',xlims');
'';



