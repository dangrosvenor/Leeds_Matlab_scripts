

%% Generic plot
nsub=nsub+1;

for i=1:length(xdat_import)
    xdat_save(i,nsub) = xdat_import(i).x;
    ydat_save(i,nsub) = ydat_import(i).y;
end


%idate_ticks_fix=0;




%run this script to get cdan(1:9) and markers(1:3) and pdan(1:3)
%LInestyles_etc




subplotting = subplotting_DRIVER;
iaxis_square=0;

for idat=1:99
    ismooth_x_import(idat)=0;
    ismooth_y_import(idat)=0;
end
idat_driver=0;



xdat_Nd = xdat_import;
%save the xdat data for use in teh RWP plot, since the Nd values don't
%align with the RWP ones

ydat_LWP = ydat_import;

%also save the marker style and colour
marker_style_overall = marker_style;
line_colour_overall = line_colour;
line_pattern_overall = line_pattern;


%% ---  Main script to do plots and save
savedir = savedir_driver;
%

graph=0; %graph choice in watervap
DRIVER_lineplot_watervap

for i=idat_micro   %length(h)-2
    set(h(i).h,'linestyle','none'); %keep this, as removes line from the legend
end

%set(h(1).h,'marker','*');
%set(h(1).h,'markersize',20);
%uistack(h(1).h,'top');

%Joint together the REMSS satellite values with a line.
x_all=[]; y_all=[];
for i=idat_micro  %1:length(xdat)-2
   x_all = cat(2,x_all,xdat(i).x);
   y_all = cat(2,y_all,ydat(i).y);   
end
[x_all,I]=sort(x_all);
y_all=y_all(I);

plot(x_all,y_all,'b','linewidth',3);

% for i=idat_micro(end)+1:length(xdat)
%     set(h(i).h,'color',line_colour(i).c);             
%     set(h(i).h,'marker',marker_style(i).m,'markerEdgeColor',line_colour(i).c,'markerFaceColor',line_colour(i).c);
%     set(h(i).h,'linestyle',line_pattern(i).p);
% 
% end






%Change the size of the window
%pos=get(gcf,'position');
%set(gcf,'position',[pos(1) pos(2) 1200 pos(4)]); 

pos=get(gca,'position');
%set(gca,'position',[pos(1) pos(2) pos(3) pos(4)*1.3]); 

savename_overall = savename;


if isave_plot==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1);
end



