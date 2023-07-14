%function [dat_PI,dat_ukesm,dat_ukesm_DJF,dat_ukesm_JJA,dat_ukesm_MAM,dat_ukesm_SON]=...
%    ACSIS_Robson_paper_load_data_generic_FUNC(var_ukesm,fscale,expt_str,expt_str2);
function [dat_PI,dat_ukesm,dat_ukesm_DJF,dat_ukesm_JJA,dat_ukesm_MAM,dat_ukesm_SON]=...
    ACSIS_Robson_paper_load_data_generic_FUNC(load_file,load_file_PI,var_ukesm,fscale);

%load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str '_all_' var_ukesm '.mat'];
%load_file_PI = load_file;

% dat_PI = load(load_file_PI,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
% dat_PI.dat_annual_ens = fscale*dat_PI.dat_annual_ens;
% dat_PI.dat_annual = fscale*dat_PI.dat_annual;
% 
% dat_ukesm = load(load_file,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
% dat_ukesm.dat_annual_ens = fscale*dat_ukesm.dat_annual_ens;
% dat_ukesm.dat_annual = fscale*dat_ukesm.dat_annual;
% 
% dat_ukesm_DJF=load(load_file,'dat_annual_ens_DJF','years_ukesm_1d','dat_annual_DJF','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
% dat_ukesm_DJF.dat_annual_ens = fscale*dat_ukesm_DJF.dat_annual_ens_DJF;
% dat_ukesm_DJF.dat_annual = dat_ukesm_DJF.dat_annual_DJF;
% 
% dat_ukesm_MAM=load(load_file,'dat_annual_ens_MAM','years_ukesm_1d','dat_annual_MAM','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
% dat_ukesm_MAM.dat_annual_ens = dat_ukesm_MAM.dat_annual_ens_MAM;
% dat_ukesm_MAM.dat_annual = dat_ukesm_MAM.dat_annual_MAM;
% 
% dat_ukesm_JJA=load(load_file,'dat_annual_ens_JJA','years_ukesm_1d','dat_annual_JJA','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
% dat_ukesm_JJA.dat_annual_ens = dat_ukesm_JJA.dat_annual_ens_JJA;
% dat_ukesm_JJA.dat_annual = dat_ukesm_JJA.dat_annual_JJA;
% 
% dat_ukesm_SON=load(load_file,'dat_annual_ens_SON','years_ukesm_1d','dat_annual_SON','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
% dat_ukesm_SON.dat_annual_ens = dat_ukesm_SON.dat_annual_ens_SON;
% dat_ukesm_SON.dat_annual = dat_ukesm_SON.dat_annual_SON;



dat_PI = load(load_file_PI,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_PI.dat_annual_ens = fscale*dat_PI.dat_annual_ens;
dat_PI.dat_annual = fscale*dat_PI.dat_annual;

dat_ukesm=load(load_file,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_ukesm.dat_annual_ens = fscale*dat_ukesm.dat_annual_ens;
dat_ukesm.dat_annual = fscale*dat_ukesm.dat_annual;
dat_ukesm.fscale = fscale;
dat_ukesm.gcm_area_UM = calc_area_lat_lon2d(dat_ukesm.gcm_Plat2D_edges_UM,dat_ukesm.gcm_Plon2D_edges_UM);

dat_ukesm_DJF=load(load_file,'dat_annual_ens_DJF','years_ukesm_1d','dat_annual_DJF','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_ukesm_DJF.dat_annual_ens = fscale*dat_ukesm_DJF.dat_annual_ens_DJF;
dat_ukesm_DJF.dat_annual = fscale*dat_ukesm_DJF.dat_annual_DJF;
dat_ukesm_DJF.gcm_area_UM = calc_area_lat_lon2d(dat_ukesm_DJF.gcm_Plat2D_edges_UM,dat_ukesm_DJF.gcm_Plon2D_edges_UM);

dat_ukesm_MAM=load(load_file,'dat_annual_ens_MAM','years_ukesm_1d','dat_annual_MAM','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_ukesm_MAM.dat_annual_ens = fscale*dat_ukesm_MAM.dat_annual_ens_MAM;
dat_ukesm_MAM.dat_annual = fscale*dat_ukesm_MAM.dat_annual_MAM;
dat_ukesm_MAM.gcm_area_UM = calc_area_lat_lon2d(dat_ukesm_MAM.gcm_Plat2D_edges_UM,dat_ukesm_MAM.gcm_Plon2D_edges_UM);

dat_ukesm_JJA=load(load_file,'dat_annual_ens_JJA','years_ukesm_1d','dat_annual_JJA','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_ukesm_JJA.dat_annual_ens = fscale*dat_ukesm_JJA.dat_annual_ens_JJA;
dat_ukesm_JJA.dat_annual = fscale*dat_ukesm_JJA.dat_annual_JJA;
dat_ukesm_JJA.gcm_area_UM = calc_area_lat_lon2d(dat_ukesm_JJA.gcm_Plat2D_edges_UM,dat_ukesm_JJA.gcm_Plon2D_edges_UM);

dat_ukesm_SON=load(load_file,'dat_annual_ens_SON','years_ukesm_1d','dat_annual_SON','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_ukesm_SON.dat_annual_ens = fscale*dat_ukesm_SON.dat_annual_ens_SON;
dat_ukesm_SON.dat_annual = fscale*dat_ukesm_SON.dat_annual_SON;
dat_ukesm_SON.gcm_area_UM = calc_area_lat_lon2d(dat_ukesm_SON.gcm_Plat2D_edges_UM,dat_ukesm_SON.gcm_Plon2D_edges_UM);
