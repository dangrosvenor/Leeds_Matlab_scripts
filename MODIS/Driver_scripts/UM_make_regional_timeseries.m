function [timser] = UM_make_regional_timeseries(dat,nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM)

[iregion_lin,iregion_lin_edges,dat_regional]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat);
timser = meanNoNan(meanNoNan(dat_regional,1),1);
