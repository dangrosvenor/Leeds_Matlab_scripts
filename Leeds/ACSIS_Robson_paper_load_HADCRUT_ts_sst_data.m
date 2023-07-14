function ACSIS_Robson_paper_load_Norris_HADCRUT_ts_sst_data()

%Load in monthly surface temperature data from HadCRUT4 dataset.

%file_dir = '/home/disk/eos10/d.grosvenor/Norris_Nature_2016_dat/'; %Disk lost... not sure if it's coming back...
file_dir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/HadCRUT_surface_temps/';
file_nc = [file_dir 'HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc'];
file_nc_baseline = [file_dir 'absolute_v5.nc'];
save_file = [file_dir 'HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_regridded_UKESM1.mat'];

%% anomalies
% Using mex NetCDF reader
% nc  = netcdf(file_nc);
% time = nc{'time'}(:); %days since 1850-1-1 00:00:00
% lat= nc{'latitude'}(:); dlat=abs(mean(diff(lat)));
% lon = nc{'longitude'}(:); dlon=abs(mean(diff(lon(1:2)))); %0 to 360 convention whereas we want -180 to +180.
% anoms = nc{'tas'}(:); % temperature anomaly
% nc=close(nc);

%Using ncread since mex can't read these newer NetCDF files and can't seem
%to convert on Olympus.
time = ncread(file_nc,'time'); %days since 1850-1-1 00:00:00 - monthly from Jan 1850 to May 2022.
lat= ncread(file_nc,'latitude'); dlat=abs(mean(diff(lat)));
lon = ncread(file_nc,'longitude'); dlon=abs(mean(diff(lon(1:2)))); %0 to 360 convention whereas we want -180 to +180.
anoms = ncread(file_nc,'tas_mean'); % temperature anomaly
baseline = ncread(file_nc_baseline,'tem'); % temperature mean for baseline (12 time values - annual clims for each month)


anoms = permute(anoms,[2 1 3]);
baseline = permute(baseline,[2 1 3]);

[gcm_Plon2D_HadCRUT_360,gcm_Plat2D_HadCRUT_360] = meshgrid(lon,lat);   
i180 = find(lon>180);
lon(i180) = lon(i180)-360;
[gcm_Plon2D_HadCRUT,gcm_Plat2D_HadCRUT] = meshgrid(lon,lat);   



%inan = find(anoms<-1e29); %NaN values seem to be inserted automatically by
%ncread
%anoms(inan) = NaN;



years = 1850:2021;
nyears = length(years);
for iyear=1:nyears
    istart = (iyear-1)*12+1;
    inds = istart:istart+11;
    tas(:,:,inds) = anoms(:,:,inds) + baseline; %add the baseline on
    tas_HadCRUT_annual(:,:,iyear) = meanNoNan(tas(:,:,inds),3);
end


UM_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_SW_up_TOA.mat';
load(UM_file,'gcm_Plon2D_UM','gcm_Plat2D_UM');
i180 = find(gcm_Plon2D_UM<0);
gcm_Plon2D_UM_360 = gcm_Plon2D_UM;
gcm_Plon2D_UM_360(i180) = gcm_Plon2D_UM_360(i180) + 360;

nt = size(tas,3);
nlat = size(gcm_Plon2D_UM,1);
nlon = size(gcm_Plon2D_UM,2);
tot_tas_anom_HadCRUT = NaN * ones([nlat nlon nt]);
tot_tas_anom_HadCRUT_annual = NaN * ones([nlat nlon nyears]);

%convert the  annual averages
for itime=1:nyears
    dat = squeeze(tas_HadCRUT_annual(:,:,itime));            
    tot_tas_anom_HadCRUT_annual(:,:,itime) = griddata(gcm_Plat2D_HadCRUT,gcm_Plon2D_HadCRUT,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);        
    %tot_tas_anom_HadCRUT_annual(:,:,itime) = griddata(gcm_Plat2D_HadCRUT_360,gcm_Plon2D_HadCRUT_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);    
end

%convert all the data
for itime=1:nt
    dat = squeeze(tas(:,:,itime));            
    %tot_tas_anom_HadCRUT(:,:,itime) = griddata(gcm_Plat2D_HadCRUT_360,gcm_Plon2D_HadCRUT_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);    
    tot_tas_anom_HadCRUT(:,:,itime) = griddata(gcm_Plat2D_HadCRUT,gcm_Plon2D_HadCRUT,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);  
end


%% save all vars in the save file
save(save_file,'-V7.3');

