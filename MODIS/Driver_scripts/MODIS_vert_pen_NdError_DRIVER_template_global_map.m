%Stripped down template for running plot_global_maps to produce a global
%map

% Set some stuff

%Define lat and lon grids
filename = '/home/disk/eos8/d.grosvenor/mat_files_various/SAVED_ann2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat';
load(filename,'MLAT','MLON');

[gcm_Plon2D_UM,gcm_Plat2D_UM] = meshgrid(MLON,MLAT);
[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
daynum_timeseries3_UM = 1; %[1:length(time)];
gcm_time_UTC_UM = 1; %[1:length(time)];
gcm_str_DRIVER='UM';  %N.B. - have to set it to this and not another name - doesn't work if do (code would need changing)

%Data to plot
%dat_modis = reshape([1:length(MLAT)*length(MLON)],[length(MLAT) length(MLON)]);
%titlenam_driver = ['Number of days with data'];
units_str_plot = '#';

LAT_val_DRIVER = [-1e9 1e9];
LON_val_DRIVER = [-1e9 1e9];

%LAT_val_DRIVER = [-10 20];
%LON_val_DRIVER = [-50 0];

irestrict_domain_DRIVER=0;

%--- run the file to set up the defaults
        plot_global_maps_defaults
        
    savedir =  savedir_DRIVER;   
        

iset_min_clim=0;
clim_min=0;
iset_max_clim=0;
clim_max=64800;

isave_plot=0;
iplot_markers=0;        
        
irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not

thresh_LAT = LAT_val_DRIVER;
thresh_LON = LON_val_DRIVER;

%--- set some options for these particular plot loops
set_screening = {'none'};
modis_data_plot = 'Map of 2D data from outside driver script';





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

if subplotting==1
    inew_figure=0;
end
                
%---  Run plot script and save
plot_global_maps
%-------------------------------

if subplotting==1
    %increase_font_size_map_figures
end

if isave_plot==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0);
    close(gcf);
end


%% Also calculate the mean values over the different regions.
ilat=0;
for iregion=1:length(iregions)    
    ireg=find(i2==iregion);
    dat_reg=[];
    dat_std_reg=[];
    dat_N_reg=[];
    dat_me_reg=[];
    
    %Loop over multiple boxes that may constitute a region (e.g. Alaska
    %region)
    for i=1:length(ireg)                
        ilat = ilat + 1;          
        LAT_val_DRIVER = LATs{ilat}; LON_val_DRIVER = LONs{ilat};
    
        ilat2 = find(MLAT >= LAT_val_DRIVER(1) & MLAT <= LAT_val_DRIVER(2));
        ilon2 = find(MLON >= LON_val_DRIVER(1) & MLON <= LON_val_DRIVER(2));  
        
        dat = dat_modis(ilat2,ilon2);        
        dat_reg = cat(1,dat_reg,dat(:));
        if exist('stdev')
            %stdev, me and nnums calculated in calling script
            dat_std = stdev(ilat2,ilon2);
            dat_std_reg = cat(1,dat_std_reg,dat_std(:)); %This will constain all of the std-devs over time for each column
            dat_N = nnums(ilat2,ilon2);
            dat_N_reg = cat(1,dat_N_reg,dat_N(:));            
            dat_me = me(ilat2,ilon2);
            dat_me_reg = cat(1,dat_me_reg,dat_me(:));                        
        else
            dat_std_reg = NaN;
            dat_N_reg = NaN;
            dat_me_reg = NaN;
        end
        
    end
    
    mean_reg_val(iregion) = meanNoNan(dat_reg,1); %mean of dat_modis data - for std plots this is
    %the mean over all locations of the relative std devs of errors (taken over time)
    max_reg_val(iregion) = maxALL(dat_reg); %mean of dat_modis data - for std plots this is    
    min_reg_val(iregion) = minALL(dat_reg); %mean of dat_modis data - for std plots this is        
    
    mean_abs_err_reg_val(iregion) = meanNoNan(dat_me_reg,1); %
    mean_std_equal_weighted_reg_val(iregion) = meanNoNan(dat_std_reg,1); %Each column is equally weighted for the mean of the std. devs
    
    %The above count each region equally no matter how many datapoints go
    %into the timeseries
    
    %This weights each column's result by the number of points in their
    %timeseries - whereas the mean_reg_val does a straight spatial mean of
    %each time-averaged value.
    [std_reg_val(iregion), me_reg_val(iregion)]=std_combine2(dat_std_reg',dat_me_reg',dat_N_reg');
end







% General instructions for case 'Generic plot specified outside of script'
%Run plot_global_maps_defaults (might want to change these)
            %Then need to set e.g. P, gcm_Plat2D_CERES, gcm_Plon2D_CERES,
            % gcm_Plat2D_edges_CERES, gcm_Plon2D_edges_CERES],
            % gcm_str, daynum_timeseries3_CERES = 1; gcm_time_UTC_CERES = 1; month_amsre=1; year_amsre=1;
            % mod_data_type='AMSRE';
            % Can set these times to one, or do properly depending on what
            % time screening you want, etc.