dat_ukesm = load(load_file,'gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
[gcm_area_UM] = calc_area_lat_lon2d(dat_ukesm.gcm_Plat2D_edges_UM,dat_ukesm.gcm_Plon2D_edges_UM);
gcm_Plon2D_UM = dat_ukesm.gcm_Plon2D_UM; gcm_Plat2D_UM = dat_ukesm.gcm_Plat2D_UM;  gcm_Plon2D_edges_UM = dat_ukesm.gcm_Plon2D_edges_UM;  gcm_Plat2D_edges_UM = dat_ukesm.gcm_Plat2D_edges_UM; 
dat_PI=[]; no_PI=1;
dat_ukesm.years_ukesm_1d = years_sw_calc';


