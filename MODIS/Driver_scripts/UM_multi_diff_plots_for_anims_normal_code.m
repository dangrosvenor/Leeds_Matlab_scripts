
if isubplot_multi_anim==1
    scrsz=get(0,'ScreenSize');
    posit=[9 50 scrsz(3)/1.46 scrsz(4)/1.46];
    hf_multi=figure('position',posit);
end

for inorm_case=1:length(UM_time_out_mvar)

    if isubplot_multi_anim==1
        hf=subplot(2,1,inorm_case);
    end

    dat_modis = UM_time_out_mvar{inorm_case}.datUM{1};


    if icoarsen==1
        %reset these in case  have already smoothed once
        gcm_Plat2D_edges_UM = UM_time_out_mvar{inorm_case}.gcm_Plat2D_edges_UM;
        gcm_Plon2D_edges_UM = UM_time_out_mvar{inorm_case}.gcm_Plon2D_edges_UM;
        gcm_Plat2D_UM = UM_time_out_mvar{inorm_case}.gcm_Plat2D_UM;
        gcm_Plon2D_UM = UM_time_out_mvar{inorm_case}.gcm_Plon2D_UM;

        coarsen_method = '2d smoothing filter'; N = 128; M=N;  %Smoothing
        coarsen_method = 'gridpoint blocks';N = 64; M=N; %Coarse graining to NxM gridpoint blocks (rather than smoothing)
        %coarsen_method = 'lat lon blocks'; dlat_target=0.25; dlon_target=0.25;
        N = 16; M=N; %was 64

        coarse_grain_UM; %script for this (operates on dat_modis and the lat/lon arrays)
    end


    time_str=datestr(UM_time_in.time_specific);
%    z_str='at ';
    z_str='';
%    diff_str{inorm_case} = [UM_time_out_mvar{1}.UM_case_out.labs_UM(inorm_case).l];
%    units_str_plot = UM_map_mvar(inorm_case).var_units_str;
    units_str_plot = UM_time_out_mvar{inorm_case}.var_units_str;
    titlenam_driver = [UM_time_out_mvar{inorm_case}.var_nice_name ' (' units_str_plot ') for ' time_str z_str];



    %% Plot the plot
    % See UM_plot_global_min_code for the min code needed for this

    %--- run the file to set up the defaults
    plot_global_maps_defaults

    if isubplot_multi_anim==1
        iover_ride_plot_global=1;
        inew_figure=0;
    end

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
    %        time_str = UM_time_out_mvar.time_str;
    %        if iz>-1
    %            z_str = ['at ' num2str(z_um(iz)/1000,'%.2f') ' km '];
    %        else
    %            z_str ='';
    %        end
    %        z_str = UM_time_out_mvar.z_str;



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

    iset_min_clim=UM_map_mvar(inorm_case).iset_min_clim;  %flags for setting the colorscale limits or not
    iset_max_clim=UM_map_mvar(inorm_case).iset_max_clim;
    clim_min = UM_map_mvar(inorm_case).clim_min;
    clim_max = UM_map_mvar(inorm_case).clim_max;

    %---  Run plot script and save
    plot_global_maps
    
    if isubplot_multi_anim==1
        set(hf_multi,'name',figlab);
    end
    

    dat_modis_save{inorm_case} = dat_modis;
    
    p=get(gca,'position');
    p(4)=p(4)*1.11;
    set(gca,'position',p);

end

if UM_map.isave_plot==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0);
    close(gcf);
end



