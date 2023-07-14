%% LWP plot (weighted by SW timeseries)

%--- run the file to set up the defaults
watervap_defaults


%var_str = 'SW_up_TOA';
titlenam = '';
xlab='Droplet concentration (cm^{-3})';
ylab='LWP_{ic} (g m^{-2})';

%legacy code, allowing to choose certain runs if want to leave some out
%(e.g. leave out obs)
%inds_plot = [1:length(fileUM)];
%idat_micro = [1:length(inds_plot)]; %for joining the lines together

clear xdat_import ydat_import
idat_driver=0;
for i=inds_plot
    idat_driver=idat_driver+1;
    
    Boutle_overall_LWP_ACI_v2_get_weights
    
    %LWP one longer than SW for some reason??
    lwp_weighted  = sum( LWP_incloud_20gmsq(i).timeseries_UM(1:end-1) .* w ) ./ sum(w);
    ydat_import(idat_driver).y = lwp_weighted; %one value per line (different linestyles)
    xdat_import(idat_driver).x = meanNoNan(Nd(i).timeseries_UM(:),1);  %Nd_overall_mean(i); %one value per line (different linestyles)    
    labs_import(idat_driver) = labs_UM(i); %one value per line (different linestyles)     
    marker_style(idat_driver) = marker_styleUM(i);
    line_colour(idat_driver) = line_colourUM(i);    
    line_pattern(idat_driver) = line_patternUM(i);      
end



%% ---  Set things specific for this script
ichoose_styles=1;
nmark=-1;
marksize=15;
lor=-1; %-99; %-99=no lengend, 1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
izlim=1;
zmin=20; zmax=130;
xlims=1;
xlimits=[4 1e3];
isave_plot=0; %This will be overwritten later
x_axis_type = 'log10_matlab';



%plotting function commands :-
Boutle_overall_LWP_ACI_v2_make_subplot_v2



