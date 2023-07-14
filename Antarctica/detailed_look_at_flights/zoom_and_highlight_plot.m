

if exist('i3D') & iplot_3D==1
    if i3D==1
        axis([xlims3D ylims3D]);
    end
    clear i3D
else
%    set(gca,'xlim',xlimits);
    if ix_distance==0
        datetick('x',datetick_type,'keepticks'); %redraws the ticks
    end
end

if iadd_terrain2==1
    add_terrain
end

if iplot_latlon2==1
    plot_latlon_lines2;
end

if length(time_highlight_path)>0
switch highlight_type

    case 'rectangle'
        x=time_highlight_path(1);

        y=ylims(1);
        h=ylims(2)-ylims(1);
        w=time_highlight_path(2)-time_highlight_path(1);

        X=[x x+w x+w x];
        Y=[y y y+h y+h];

        %        rectangle('Position',[x/24,y,w/24,h],'LineWidth',2,'LineStyle','--');
        line(X([1 4])/24,Y([1 4]),'LineWidth',2,'LineStyle','--','color','k');
        line(X([2 3])/24,Y([2 3]),'LineWidth',2,'LineStyle','--','color','k');
end
end

    hold on

if isave_highlight==1    
    iremove=findstr(savename,':');
    savename(iremove)='_';
    
    iremove=findstr(savename,'\');
    savename(iremove)='';    
    
    iremove=findstr(savename,'/');
    savename(iremove)='';
    
    iremove=findstr(savename,'>');
    savename(iremove)='';
    
    iremove=findstr(savename,'<');
    savename(iremove)='';
    
    print(gcf,[savedir savename '.emf'],'-dmeta');    
%    print(gcf,[savedir savename '.eps'],'-depsc');    
    print(gcf,[savedir savename '.tiff'],'-dtiff');
    saveas(gcf,[savedir savename],'fig');
    
    close(gcf);
end



