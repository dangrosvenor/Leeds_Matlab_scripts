file01='/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/r19i1p1f2/Nd_clw_weighted_ESGF/cdnc_AERmon_UKESM1-0-LL_historical_r19i1p1f2_gn_185001-189912.nc_197001011500_Nd_clw_weighted_ESGF_total_column_to_zdomain_top_saved.nc3';
nc=netcdf(file01);
dat = nc{'Nd_clw_weighted_ESGF'}(1,:,:);
time = nc{'time'}(:);
nc=close(nc);


file01='/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/r19i1p1f2/Nd_clw_weighted_ESGF/cdnc_AERmon_UKESM1-0-LL_historical_r19i1p1f2_gn_190001-194912.nc_197201211500_Nd_clw_weighted_ESGF_total_column_to_zdomain_top_saved.nc3';
nc=netcdf(file01);
%dat = nc{'Nd_clw_weighted_ESGF'}(1,:,:);
time02 = nc{'time'}(:);
nc=close(nc);

file01='/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/r19i1p1f2/Nd_clw_weighted_ESGF_no_dz/cdnc_AERmon_UKESM1-0-LL_historical_r19i1p1f2_gn_185001-189912.nc_197001011500_Nd_clw_weighted_ESGF_no_dz_total_column_to_zdomain_top_saved.nc3';
nc=netcdf(file01);
dat_no_dz = nc{'Nd_clw_weighted_ESGF_no_dz'}(1,:,:);
time_no_dz = nc{'time'}(:);
nc=close(nc);

file01='/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/r19i1p1f2/Nd_clw_weighted_ESGF/output_pre_mask_etc_changes/cdnc_AERmon_UKESM1-0-LL_historical_r19i1p1f2_gn_185001-189912.nc_197001011500_Nd_clw_weighted_ESGF_total_column_to_zdomain_top_saved.nc3';
nc=netcdf(file01);
dat_old = nc{'Nd_clw_weighted_ESGF'}(1,:,:);
time_old = nc{'time'}(:);
nc=close(nc);

%%

clear dat_cnrm_cm6
for i=1:3

file01=['/home/disk/eos15/d.grosvenor/ESGF/CNRM-CERFACS/CNRM-CM6-1/r' istr 'i1p1f2/output/Nd_clw_weighted_ESGF_no_dz/cdnc_AERmon_CNRM-CM6-1_historical_r' istr 'i1p1f2_gr_185001-201412.nc_197001011530_Nd_clw_weighted_ESGF_no_dz_total_column_to_zdomain_top_saved.nc3'];
nc=netcdf(file01);
dat_cnrm_cm6(:,:,i) = nc{'Nd_clw_weighted_ESGF_no_dz'}(1,:,:);
time_cnrm = nc{'time'}(:);
nc=close(nc);

end

dat_cnrm_cm6_av = meanNoNan(dat_cnrm_cm6,3);

qpcolor(dat_cnrm_cm6_av/1e6);
title('CNRM-CM6-1');


%%
clear  dat_cnrm_esm
for i=1:3
    istr = num2str(i);
    file01=['/home/disk/eos15/d.grosvenor/ESGF/CNRM-CERFACS/CNRM-ESM2-1/historical/r' istr 'i1p1f2/output/Nd_clw_weighted_ESGF_no_dz/cdnc_AERmon_CNRM-ESM2-1_historical_r' istr 'i1p1f2_gr_185001-201412.nc_197001011530_Nd_clw_weighted_ESGF_no_dz_total_column_to_zdomain_top_saved.nc3'];
    nc=netcdf(file01);
    dat_cnrm_esm(:,:,i) = nc{'Nd_clw_weighted_ESGF_no_dz'}(1,:,:);
    time_cnrm = nc{'time'}(:);
    nc=close(nc);
end

dat_cnrm_esm_av = meanNoNan(dat_cnrm_esm,3);

qpcolor(dat_cnrm_esm_av/1e6);
title('CNRM-ESM2-1');



