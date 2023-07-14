try

iUTC=1;

gcm_str_DRIVER='UM';  %N.B. - have to set it to this and not another name - doesn't work if do (code would need changing)
%[gcm_Plon2D_UM,gcm_Plat2D_UM] = meshgrid(MLON,MLAT);

daynum_timeseries3_UM = 1; %[1:length(time)];
gcm_time_UTC_UM = 1; %[1:length(time)];

if ~exist('igeneric_plot')
    igeneric_plot=0;
end

if icoarse_grain==1
    gcm_str_DRIVER='GENERIC';        
        
    dat_modis = reduce_matrix_subsample_mean(dat_modis,M_coarse_grain,N_coarse_grain);
    gcm_Plat2D_GENERIC = reduce_matrix_subsample_mean(gcm_Plat2D_UM,M_coarse_grain,N_coarse_grain);
    gcm_Plon2D_GENERIC = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse_grain,N_coarse_grain);
    [gcm_Plat2D_edges_GENERIC,gcm_Plon2D_edges_GENERIC] = get_edges_lat_lon(gcm_Plat2D_GENERIC,gcm_Plon2D_GENERIC);
    daynum_timeseries3_GENERIC = daynum_timeseries3_UM;
    gcm_time_UTC_GENERIC = gcm_time_UTC_UM;
    
elseif igeneric_plot==1
    gcm_str_DRIVER='GENERIC'; 
end


set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.

isquare=0;

% if iUTC==1
%     time_shift=0;
%     time_format_str=' UTC';
% else
%     time_format_str=' LST';
% end




clear clim_min clim_max

%Script for variable labels and clims, units etc.
UM_var_defs
if ~exist('titlenam_driver')
    titlenam_driver = remove_character(var_UM,'_',' ');
end


units_str_plot = '';

%%

if ~exist('ioverride_LAT_plots') | ioverride_LAT_plots==0

LAT_val_DRIVER = [-1e9 1e9]; LON_val_DRIVER = [-1e9 1e9]; %Can restrict the plot region with this if set irestrict_domain_DRIVER=1

LAT_val_DRIVER = [-40 0]; LON_val_DRIVER = [-100 -60]; 

% LAT_val_DRIVER = [-24 -16]; LON_val_DRIVER = [-80 -72]; 

LAT_val_DRIVER = [-30 -10]; LON_val_DRIVER = [-86 -66]; %Close up of region surrounding the UM domain

LAT_val_DRIVER = [20 80]; LON_val_DRIVER = [-86 0]; %N Atlantic

LAT_val_DRIVER = [0 80]; LON_val_DRIVER = [-179 179]; %N Atlantic
LAT_val_DRIVER = [-10 60]; LON_val_DRIVER = [-75 0]; %N Atlantic

LAT_val_DRIVER = [-10 60]; LON_val_DRIVER = [-85 0]; %N Atlantic2
%LAT_val_DRIVER = [27 62]; LON_val_DRIVER = [-72 0]; %N Atlantic animation global and nested
LAT_val_DRIVER = [-10 80]; LON_val_DRIVER = [-85 0]; %N Atlantic - looking further north
LAT_val_DRIVER = [-10 80]; LON_val_DRIVER = [-85 30]; %N Atlantic - looking further north

%LAT_val_DRIVER = [30 45]; LON_val_DRIVER = [-60 -10]; %N Atlantic region of high negative forcing ocean only

%LAT_val_DRIVER = [-30 30]; LON_val_DRIVER = [-43 43]; %Africa

%LAT_val_DRIVER = [-30 0]; LON_val_DRIVER = [-100 -60]; %VOCALS
%LAT_val_DRIVER = [-24 -18]; LON_val_DRIVER = [-80 -72]; %VOCALS

%LAT_val_DRIVER = [-40 0]; LON_val_DRIVER = [-176 176]; %Rosenfeld ALL
%LAT_val_DRIVER = [-40 0]; LON_val_DRIVER = [-176 -40]; %Rosenfeld to equator

%LAT_val_DRIVER = [40 82]; LON_val_DRIVER = [-85 55]; %Iceland N Atlantic region

%LAT_val_DRIVER = [-90 90]; LON_val_DRIVER = [-160 160];
%LAT_val_DRIVER = [-30 0]; LON_val_DRIVER = [-100 -60];

%LAT_val_DRIVER = [-80 0]; LON_val_DRIVER = [60 179]; %SOCRATES region.

else
    
    LAT_val_DRIVER = LAT_val_DRIVER_override;
    LON_val_DRIVER = LON_val_DRIVER_override;

end

%irestrict_domain_DRIVER=1;

%--- run the file to set up the defaults
        plot_global_maps_defaults
        
        icontour = icontour_DRIVER;
        cont_col_str = cont_col_str_DRIVER;
        iplot_mgrid_lines = iplot_mgrid_lines_DRIVER;
        ioverride_ticks = ioverride_ticks_DRIVER;
        if exist('ioverride_proj_type') & ioverride_proj_type==1
           proj_type = proj_type_DRIVER;                    
        end
        
irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not

thresh_LAT = LAT_val_DRIVER;
thresh_LON = LON_val_DRIVER;

%--- set some options for these particular plot loops
set_screening = {'none'};
modis_data_plot = 'Map of 2D data from outside driver script';

if exist('clim_min')
iset_min_clim=1;
%clim_min=0;
iset_max_clim=1;
%clim_max=64800;
end

%isave_plot=0;
iplot_markers=0;



mod_data_type='AMSRE';
gcm_str_select = gcm_str_DRIVER;
%        daynum_timeseries3_UM = [1:length(time)];
%        gcm_time_UTC_UM = [1:length(time)];


ifull_swath=0;
igcm_screen=0;



inew_figure = 0;


%--- Apply override flags
ioverride_plotglobal_thresh=1; %Override most of the options (what to plot, etc.)
% iocean_only=1;
ioverride_time_selection=0; %Override the times to include
ioverride_plotglobal_loc=1; %Override the location of the plot window
ioverride_years_time_screen=0; %Override years for screening?
iover_ride_plot_global=1; %overrides inew_figure=1; supress_colorbar=0; i_increase_font_size_map_figures_OFF = 0;
%(all set in plot_global_maps_defaults)


%---  Run plot script and save
plot_global_maps
%-------------------------------

clear ilabel_colorbar

lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));

%tit_str_clean = [titlenam_driver];
% tit_str = [titlenam_driver ' mean=' num2str(Pmean,'%.2g')];
% %tit_str='test';
% ncols=55; nlines=4; %nlines is the no. of lines we are aiming for (will add blank ones to make up to this later)
% tit_wrapped = textwrap({tit_str},ncols);
% %Add extra lines if they don't exist
% for iw=length(tit_wrapped)+1:nlines
%     tit_wrapped{iw}='';
% end
% tit_wrapped{iw} = tit_str_clean;
% title(tit_wrapped);

add_str = var_UM;
clear title_full       
    title_full{1}=[titlenam_driver add_str ', mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str ];
    title(title_full);

tit_str_clean='';




savename=['/home/disk/eos1/d.grosvenor/modis_work/lwp_vs_Nd_plots/global_' titlenam_driver];

if isave_plot==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_LWP_' titlenam_driver '_vs_satellite'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%    close(gcf);
end



% General instructions for case 'Generic plot specified outside of script'
%Run plot_global_maps_defaults (might want to change these)
            %Then need to set e.g. P, gcm_Plat2D_CERES, gcm_Plon2D_CERES,
            % gcm_Plat2D_edges_CERES, gcm_Plon2D_edges_CERES],
            % gcm_str, daynum_timeseries3_CERES = 1; gcm_time_UTC_CERES = 1; month_amsre=1; year_amsre=1;
            % mod_data_type='AMSRE';
            % Can set these times to one, or do properly depending on what
            % time screening you want, etc.
            
            
            
if iplot_wind_arrows==1            
            
            
nx_quiver=25; %number of arrows to draw for x
ny_quiver=25; %number for y
scale_speed_quiver = [15 15]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
%scale_speed_quiver = [7 7]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
%scale_speed_quiver = [10 10]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
scale_speed_quiver = [1 1]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
scale_speed_quiver = [0.5 0.5]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed



nx_quiver=100; %number of arrows to draw for x
ny_quiver=100; %number for y

nx_quiver = size(gcm_Plon2D_UM,2)/4; %All of th arrows divided by a factor
ny_quiver = size(gcm_Plon2D_UM,1)/4;

%nx_quiver=500; %number of arrows to draw for x
%ny_quiver=500; %number for y

%nx_quiver=25; %number of arrows to draw for x
%ny_quiver=25; %number for y      


%[ax_orig]=plot_wind_quiver_arrows(u_quiver,v_quiver,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');

%m_quiver(280,-20,U_glm(300,300),V_glm(300,300),0,'w');

%m_quiver(gcm_Plon2D_UM,gcm_Plat2D_UM,U_glm,V_glm(1:end-1,:));



%convert the lat and lon values to map coords
[x_quiver,y_quiver]=m_ll2xy(gcm_Plon2D_UM,gcm_Plat2D_UM);
%[ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');

[ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');    
%[ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,'k');    


%Plot the UM nest domain
ioverride_box_colour=1;
box_lwidth = 3;
LAT_val = [-22.70 -17.28]; LON_val =[-78.93 -73.08]; %12th Nov UM FULL domain
col_str='k-';
itext_in_box=0;

end

%plot_box_on_map

%clear ioverride_LAT_plots
clear ioverride_proj_type
catch err_catch
%    clear ioverride_LAT_plots
    clear ioverride_proj_type
    rethrow(err_catch)
end