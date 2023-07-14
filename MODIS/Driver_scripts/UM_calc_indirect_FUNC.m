function UM_calc_indirect_FUNC(SW_var,SW_surf_or_TOA


% For indirect effect we do clean minus clear_clean.
% I.e., no aerosol effect with cloud effect - no aerosol effect without cloud effect
    %no direct or indirect effects
    
 var_UM = [SW_var '_clean_clear_'  SW_surf_or_TOA]; 
    
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [SW_clean_clear_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [SW_clean_clear_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    
    cloud_effect_PI_ALL = SW_nodirect_PI_ALL - SW_clean_clear_PI_ALL;    
    cloud_effect_PD_ALL = SW_nodirect_PD_ALL - SW_clean_clear_PD_ALL;
    indirect_ALL = cloud_effect_PD_ALL - cloud_effect_PI_ALL;