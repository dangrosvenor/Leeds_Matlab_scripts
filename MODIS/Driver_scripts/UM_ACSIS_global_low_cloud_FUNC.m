function [low_CF_PD_ALL,nT] = UM_ACSIS_global_low_cloud_FUNC(cloud_height_str,cloud_input,um_case_PD,load_type,time_inds,gcm_Plat2D_UM)

switch cloud_input
    case 'UM'
        var_UM = [cloud_height_str '_cloud_amount'];
        icosp_mask=0;
    case 'CALIPSO'
        var_UM = ['calipso_' cloud_height_str '_cloud_amount'];
        var_UM_mask = [var_UM '_mask'];
        icosp_mask=1;
    case 'MODIS'
        var_UM = ['modis_' cloud_height_str '_cloud_amount'];
        var_UM_mask = ['modis_misr_issp_cloud_amount_mask'];
        icosp_mask=1;
end
 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [low_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    low_CF_PD_ALL(low_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??  
    
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [low_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        low_CF_PD_ALL(low_CF_mask==0) = NaN;
    end
    