%% Plot UKESM ensemble mean values for May for comparison with Leighton (used May 2015).

irestrict_domain_DRIVER=0;

clear time_choice; time_choice.time_range = [datenum('01-May-2009') datenum('31-May-2009')];
[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

var_UM = 'high_cloud_amount';
um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
[high_CF,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
high_CF(high_CF>1.001)=NaN; %got some strange spikes in the data for mid-cloud??

dat_modis = meanNoNan(high_CF,3); var_UM = 'High CF May 2009';
    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 1]);
