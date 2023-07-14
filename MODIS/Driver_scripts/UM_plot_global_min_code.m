


        %--- run the file to set up the defaults
        plot_global_maps_defaults

        irestrict_domain=UM_time_in.irestrict_domain_DRIVER; %whether to restrict the domain or not


        thresh_LAT = UM_time_in.LAT_val_DRIVER;
        thresh_LON = UM_time_in.LON_val_DRIVER;

        %--- set some options for these particular plot loops
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';

        %         iset_min_clim=1;
        %         clim_min=0;
        %         iset_max_clim=1;
        %         clim_max=300;
        %
        %         isave_plot=0;
        %         iplot_markers=0;



        %Set various things
        %        time = time_matlab(it_driver);
        %Round to the nearest minute as sometimes get 18:59:59
        %        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM');
        %        time_str = UM_time_out.time_str;
        %        if iz>-1
        %            z_str = ['at ' num2str(z_um(iz)/1000,'%.2f') ' km '];
        %        else
        %            z_str ='';
        %        end
        %        z_str = UM_time_out.z_str;



        mod_data_type='AMSRE';
        gcm_str_select='UM';
        %        daynum_timeseries3_UM = [1:length(time)];
        %        gcm_time_UTC_UM = [1:length(time)];

        %        gcm_Plon2D_UM = lons2D_driver;
        %        gcm_Plat2D_UM = lats2D_driver;

        %        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(lats2D_driver,lons2D_driver);



        %        i_dpcolor=1;
        ifull_swath=0;
        igcm_screen=0;
        gcm_time_matlab_UM = 1;
        gcm_time_UTC_UM = 1;
        daynum_timeseries3_UM = 1;
        modisyear_timeseries3_UM = 1;




        %--- Apply override flags
        ioverride_plotglobal_thresh=1; %Override most of the options (what to plot, etc.)
        % iocean_only=1;
        ioverride_time_selection=0; %Override the times to include
        ioverride_plotglobal_loc=1; %Override the location of the plot window
        ioverride_years_time_screen=0; %Override years for screening?

        iset_min_clim=UM_map.iset_min_clim;  %flags for setting the colorscale limits or not
        iset_max_clim=UM_map.iset_max_clim;
        clim_min = UM_map.clim_min;
        clim_max = UM_map.clim_max;

        %---  Run plot script and save
        plot_global_maps

