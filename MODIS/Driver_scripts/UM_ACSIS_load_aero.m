function [dat_aero_map] = UM_ACSIS_load_aero(um_case,dirUM,var_UM,run_type,time_inds,load_type,gcm_Plat2D_UM,var_UM2)


    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,[],[],var_UM2);   
    [dat_aero,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
    dat_aero_map = meanNoNan(dat_aero,3);
    