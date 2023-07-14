%Choose runs, etc. from the UM_ACSIS_SW_vs_cloud_properties_global.m script
%first.

%% scavenging vs no scav, etc., plots    

if exist('um_case_CONV_SCAV_ON')
    um_case=um_case_CONV_SCAV_ON; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM_Nd,run_type,load_type);   
    [Nd_CONV_SCAV_ON_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
    Nd_CONV_SCAV_ON_ALL_map = meanNoNan(Nd_CONV_SCAV_ON_ALL,3);
 
    
    dat_modis = (Nd_PD_map - Nd_CONV_SCAV_ON_ALL_map)/1e6; var_UM = 'Effect of large scale scav (All scav ON minus Large-scale OFF ); cm^{-3})';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]); 
    
    dat_modis = (Nd_PD_map - Nd_PI_map)/1e6; var_UM = 'Effect of convective+large-scale scav (All scav ON minus ALL scav OFF); cm^{-3})';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]); 
    
    dat_modis = ( Nd_CONV_SCAV_ON_ALL_map - Nd_PI_map )/1e6; var_UM = 'Effect of convective scav (Large-scale OFF, Conv ON minus All scav OFF minus ); cm^{-3})';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]); 
    
    um_case=um_case_r37pt5; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM_Nd,run_type,load_type);   
    [Nd_r37pt5_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
    Nd_r37pt5_map = meanNoNan(Nd_r37pt5_ALL,3);
    
    um_case=um_case_r37pt5_rain_frac1; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM_Nd,run_type,load_type);   
    [Nd_r37pt5_rain_frac1_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
    Nd_r37pt5_rain_frac1_map = meanNoNan(Nd_r37pt5_rain_frac1_ALL,3);    
    
    dat_modis = (Nd_r37pt5_map - Nd_PD_map )/1e6; var_UM = '(r=37.5 minus All scav ON); cm^{-3})';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]); 
    
    dat_modis = (Nd_r37pt5_rain_frac1_map -  Nd_r37pt5_map)/1e6; var_UM = '(r37.5, rain_frac=1 minus r=37.5); cm^{-3})';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]); 
    
end   

%% Aerosol plots

um_case=um_case_r37pt5; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
dat_modis = UM_ACSIS_load_aero(um_case,dirUM,[var_UM_aero '_z3200'],run_type,time_inds,load_type,gcm_Plat2D_UM,var_UM_aero);
dat_modis = dat_modis / 3076.664; %Since is total integrated over 3076m - per m3
dat_r37pt5 = dat_modis;

um_case=um_case_r37pt5_rain_frac1; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1;
dat_modis = UM_ACSIS_load_aero(um_case,dirUM,[var_UM_aero '_z3200'],run_type,time_inds,load_type,gcm_Plat2D_UM,var_UM_aero);
dat_modis = dat_modis / 3076.664; %Since is total integrated over 3076m
dat_r37pt5_rain_frac1 = dat_modis;

um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1;
dat_modis = UM_ACSIS_load_aero(um_case,dirUM,[var_UM_aero '_z3200'],run_type,time_inds,load_type,gcm_Plat2D_UM,var_UM_aero);
dat_modis = dat_modis / 3076.664; %Since is total integrated over 3076m
dat_standard = dat_modis;



dat_modis = (dat_r37pt5 - dat_standard )/1e6; var_UM = [var_UM_aero ' (r=37.5 minus All scav ON; cm^{-3})'];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-50 50]);

dat_modis = (dat_r37pt5_rain_frac1 -  dat_r37pt5)/1e6; var_UM = [var_UM_aero ' (r37.5, rain_frac=1 minus r=37.5; cm^{-3})'];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-50 50]);

