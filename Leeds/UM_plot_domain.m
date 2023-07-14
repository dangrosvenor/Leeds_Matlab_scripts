%Stripped down template for running plot_global_maps to produce a global
%map

% Set some stuff

%Define lat and lon grids
%filename = '/home/disk/eos8/d.grosvenor/mat_files_various/SAVED_ann2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat';
%load(filename,'MLAT','MLON');
gcm_str_DRIVER='UM';  %N.B. - have to set it to this and not another name - doesn't work if do (code would need changing)
%[gcm_Plon2D_UM,gcm_Plat2D_UM] = meshgrid(MLON,MLAT);
[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
daynum_timeseries3_UM = 1; %[1:length(time)];
gcm_time_UTC_UM = 1; %[1:length(time)];

%Data to plot
dat_modis = ones(size(gcm_Plon2D_UM));
titlenam_driver = ['UM nested domain N Atlantic (west)'];
units_str_plot = '';

LAT_val_DRIVER = [-1e9 1e9];
LON_val_DRIVER = [-1e9 1e9]; %Can restrict the plot region with this if set irestrict_domain_DRIVER=1

LAT_val_DRIVER = [20 75];
LON_val_DRIVER = [-100 20]; %Also works to show more of the globe

irestrict_domain_DRIVER=1;

%--- run the file to set up the defaults
        plot_global_maps_defaults
        
irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not

thresh_LAT = LAT_val_DRIVER;
thresh_LON = LON_val_DRIVER;

%--- set some options for these particular plot loops
set_screening = {'none'};
modis_data_plot = 'Map of 2D data from outside driver script';

iset_min_clim=0;
clim_min=0;
iset_max_clim=0;
clim_max=64800;

isave_plot=0;
iplot_markers=0;



mod_data_type='AMSRE';
gcm_str_select = gcm_str_DRIVER;
%        daynum_timeseries3_UM = [1:length(time)];
%        gcm_time_UTC_UM = [1:length(time)];


ifull_swath=0;
igcm_screen=0;




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
if isave_plot==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0);
    close(gcf);
end



% General instructions for case 'Generic plot specified outside of script'
%Run plot_global_maps_defaults (might want to change these)
            %Then need to set e.g. P, gcm_Plat2D_CERES, gcm_Plon2D_CERES,
            % gcm_Plat2D_edges_CERES, gcm_Plon2D_edges_CERES],
            % gcm_str, daynum_timeseries3_CERES = 1; gcm_time_UTC_CERES = 1; month_amsre=1; year_amsre=1;
            % mod_data_type='AMSRE';
            % Can set these times to one, or do properly depending on what
            % time screening you want, etc.