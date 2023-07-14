%% SW plot

%--- run the file to set up the defaults
watervap_defaults

graph=0; %graph choice in watervap
var_str = 'SW_up_TOA';
titlenam = '';
xlab='Droplet concentration (cm^{-3})';
ylab='SW_{\uparrow} (W m^{-2})';

inds_plot = [1:length(fileUM)];
%idat_micro = [1:length(inds_plot)]; %for joining the lines together

clear xdat_import ydat_import
idat_driver=0;
for i=inds_plot
    idat_driver=idat_driver+1;
    
    Boutle_overall_LWP_ACI_v2_get_weights
    
    ydat_import(idat_driver).y = sum( SW_up_TOA(i).timeseries_UM(:)'.* w ) ./ sum(w);
    xdat_import(idat_driver).x = meanNoNan(Nd(i).timeseries_UM(:),1);  %Nd_overall_mean(i); %one value per line (different linestyles)    
    labs_import(idat_driver) = labs_UM(i); %one value per line (different linestyles)     
    marker_style(idat_driver) = marker_styleUM(i);
    line_colour(idat_driver) = line_colourUM(i);    
    line_pattern(idat_driver) = line_patternUM(i); 
%    SW_weighted(idat_driver) = sum( SW_up_TOA(i).timeseries_UM(:)'.* w ) ./ sum(w);
end



%% ---  Set things specific for this script
ichoose_styles=1;
nmark=-1;
marksize=15;
lor=-1; %-99=no legend, 1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
izlim=1;
zmin=0; zmax=250; %For unweighted means (all times excep frist 6 hours)
zmin=0; zmax=400; %When weighted by SW_down_TOA
%zmin=0; zmax=200;
xlims=1;
xlimits=[4 1e3];
isave_plot=0; %This will be overwritten later
x_axis_type = 'log10_matlab';



%plotting function commands :-
Boutle_overall_LWP_ACI_v2_make_subplot_v2


%% not using rest of code
return



nsub=nsub+1;

%--- set some options for this particular plot
graph=0; %graph choice in watervap



%idate_ticks_fix=0;

inds_plot = [1:length(fileUM)];
idat_micro = [1:length(inds_plot)]; %for joining the lines together


%run this script to get cdan(1:9) and markers(1:3) and pdan(1:3)
%LInestyles_etc

%--- run the file to set up the defaults
watervap_defaults

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
ichoose_styles=1;
nmark=-1;
marksize=15;
lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
izlim=1;
zmin=0; zmax=110;
xlims=1;
xlimits=[4 1e3];
isave_plot=0; %This will be overwritten later
x_axis_type = 'log10_matlab';

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
pos=get(gcf,'position');
%set(gcf,'position',[pos(1) pos(2) 1200 pos(4)]); 

savename_overall = savename;


if isave_plot==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1);
end



