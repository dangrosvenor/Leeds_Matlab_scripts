function add_points_to_plot(xpos,ypos,point_labs,msize,tsize)
%function add_points_to_plot(xpos,ypos,point_labs,msize,tsize)
%plots points on the gca graph specified by xpos(i).x, ypos(i).y and point_labs(i).lab
%recommend textsize of tsize=11
%and msize=8

ylims=get(gca,'ylim');
yspan=diff(ylims);

for i=1:length(xpos)
%    plot(xpos(i).x,ypos(i).y,marker,'markersize',size);
    plot(xpos(i).x,ypos(i).y,'rs','MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',msize);
    h=text(xpos(i).x,ypos(i).y-yspan/40,[point_labs(i).lab]);
    set(h,'fontsize',tsize,'fontweight','b','horizontalalignment','center')
    
end

