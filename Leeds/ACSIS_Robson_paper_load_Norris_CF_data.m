function ACSIS_Robson_paper_load_Norris_CF_data()

%Load in monthly cloud fraction anomaly data from Norris Nature 2016 paper

%file_dir = '/home/disk/eos10/d.grosvenor/Norris_Nature_2016_dat/'; %Disk lost... not sure if it's coming back...
file_dir = '/home/disk/eos5/d.grosvenor/modis_work/ACSIS_Nd_trends/Norris_CF_data/';
save_file = [file_dir 'Norris_cf_anom_dat.mat'];

file_nc_isccp = [file_dir 'isccp.tot_amt.corr.jan1983-dec2009.nc'];
file_nc = [file_dir 'patmosx.tot_cld.corr.jan1983-dec2009.nc'];

%% ISCCP

nc  = netcdf(file_nc_isccp);

time = nc{'decimal_year'}(:); %center point time in decimal years - i.e., first data-point is 1983.0416
%mid-Jan 1983.
%time_matlab_in(ifile) = datenum('01-Jan-1970') + time;

lat= nc{'lat'}(:); dlat=abs(mean(diff(lat)));
lon = nc{'lon'}(:); dlon=abs(mean(diff(lon(1:2)))); %0 to 360 convention whereas we want -180 to +180.
[gcm_Plon2D_ISCCP_360,gcm_Plat2D_ISCCP_360] = meshgrid(lon,lat);   
i180 = find(lon>180);
lon(i180) = lon(i180)-360;
[gcm_Plon2D_ISCCP,gcm_Plat2D_ISCCP] = meshgrid(lon,lat);   

cf_ISCCP_load = nc{'total_vis_ir_amt_monthly_anom'}(:); % anomalies
nc=close(nc);

inan = find(cf_ISCCP_load>9e36);
cf_ISCCP_load(inan) = NaN;

years = 1983:2009;
nyears = length(years);
for iyear=1:nyears
    istart = (iyear-1)*12+1;
    cf_ISCCP_annual(iyear,:,:) = meanNoNan(cf_ISCCP_load(istart:istart+11,:,:),1);
end


UM_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_SW_up_TOA.mat';
load(UM_file,'gcm_Plon2D_UM','gcm_Plat2D_UM');
i180 = find(gcm_Plon2D_UM<0);
gcm_Plon2D_UM_360 = gcm_Plon2D_UM;
gcm_Plon2D_UM_360(i180) = gcm_Plon2D_UM_360(i180) + 360;

nt = size(cf_ISCCP_load,1);
nlat = size(gcm_Plon2D_UM,1);
nlon = size(gcm_Plon2D_UM,2);
totcf_anom_ISCCP = NaN * ones([nt nlat nlon]);

for itime=1:nyears
    dat = squeeze(cf_ISCCP_annual(itime,:,:));            
    totcf_anom_ISCCP_annual(itime,:,:) = griddata(gcm_Plat2D_ISCCP,gcm_Plon2D_ISCCP,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);        
    %totcf_anom_ISCCP_annual(itime,:,:) = griddata(gcm_Plat2D_ISCCP_360,gcm_Plon2D_ISCCP_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);    
end

for itime=1:nt
    dat = squeeze(cf_ISCCP_load(itime,:,:));            
    %totcf_anom_ISCCP(itime,:,:) = griddata(gcm_Plat2D_ISCCP_360,gcm_Plon2D_ISCCP_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);    
    totcf_anom_ISCCP(itime,:,:) = griddata(gcm_Plat2D_ISCCP,gcm_Plon2D_ISCCP,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);  
end


%% PATMOS

nc  = netcdf(file_nc);

time = nc{'decimal_year'}(:); %center point time in decimal years - i.e., first data-point is 1983.0416
%mid-Jan 1983.
%time_matlab_in(ifile) = datenum('01-Jan-1970') + time;

lat= nc{'lat'}(:); dlat=abs(mean(diff(lat)));
lon = nc{'lon'}(:); dlon=abs(mean(diff(lon(1:2)))); %0 to 360 convention whereas we want -180 to +180.
[gcm_Plon2D_PATMOS_360,gcm_Plat2D_PATMOS_360] = meshgrid(lon,lat);   
i180 = find(lon>180);
lon(i180) = lon(i180)-360;
[gcm_Plon2D_PATMOS,gcm_Plat2D_PATMOS] = meshgrid(lon,lat);   

cf_PATMOS_load = nc{'aft_asc_total_cld_monthly_anom'}(:); % anomalies
nc=close(nc);

inan = find(cf_PATMOS_load>9e36);
cf_PATMOS_load(inan) = NaN;

years = 1983:2009;
nyears = length(years);
for iyear=1:nyears
    istart = (iyear-1)*12+1;
    cf_PATMOS_annual(iyear,:,:) = meanNoNan(cf_PATMOS_load(istart:istart+11,:,:),1);
end


UM_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_SW_up_TOA.mat';
load(UM_file,'gcm_Plon2D_UM','gcm_Plat2D_UM');
i180 = find(gcm_Plon2D_UM<0);
gcm_Plon2D_UM_360 = gcm_Plon2D_UM;
gcm_Plon2D_UM_360(i180) = gcm_Plon2D_UM_360(i180) + 360;

nt = size(cf_PATMOS_load,1);
nlat = size(gcm_Plon2D_UM,1);
nlon = size(gcm_Plon2D_UM,2);
totcf_anom_PATMOS = NaN * ones([nt nlat nlon]);

for itime=1:nyears
    dat = squeeze(cf_PATMOS_annual(itime,:,:));            
    %totcf_anom_PATMOS_annual(itime,:,:) = griddata(gcm_Plat2D_PATMOS_360,gcm_Plon2D_PATMOS_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);    
    totcf_anom_PATMOS_annual(itime,:,:) = griddata(gcm_Plat2D_PATMOS,gcm_Plon2D_PATMOS,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);        
end

for itime=1:nt
    dat = squeeze(cf_PATMOS_load(itime,:,:));            
    %totcf_anom_PATMOS(itime,:,:) = griddata(gcm_Plat2D_PATMOS_360,gcm_Plon2D_PATMOS_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);    
    totcf_anom_PATMOS(itime,:,:) = griddata(gcm_Plat2D_PATMOS,gcm_Plon2D_PATMOS,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);
end



save(save_file,'-V7.3');

