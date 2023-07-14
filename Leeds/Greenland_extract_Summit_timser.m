
model_str = 'u-bz785, wind nudging, abv. BL'; load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_Nudged_u-bz785_all_LWP_2-391.mat';
dat_ALL=load(load_file,'dat_ens');
t = load(load_file,'years_ukesm','months_ukesm');
years = t.years_ukesm'; years = years(:); %all the years one after the other for monthly data.
months = t.months_ukesm'; months = months(:); %all the years one after the other for monthly data.
time_decimal_year = years + (months-0.5)/12;

%Summit is at 72 34' 46.50'' N 38 27' 33.07'' W according to Wikipedia
grid_dat = load(load_file,'gcm_Plat2D_UM','gcm_Plon2D_UM');
lat_summit = 72 + 34/60 + 46.5/3600;
lon_summit = - (38 + 27/60 + 33.07/3600);
[minval ilat] = min(abs(grid_dat.gcm_Plat2D_UM(:,1) - lat_summit));
[minval ilon] = min(abs(grid_dat.gcm_Plon2D_UM(1,:) - lon_summit)); %lon is in -180 to +180 form

lwp_summit = dat_ALL.dat_ens(1,:,ilat,ilon)*1e3; %g/m2
n=1;
lwp_summit_3_by_3 = dat_ALL.dat_ens(1,:,ilat-n:ilat+n,ilon-n:ilon+n)*1e3; %g/m2
lwp_summit_3_by_3_av = meanNoNan(lwp_summit_3_by_3(:,:,:),3);


model_str = 'u-bz785, wind nudging, abv. BL'; load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_Nudged_u-bz785_all_IWP_2-392.mat';
dat_ALL=load(load_file,'dat_ens');
iwp_summit = dat_ALL.dat_ens(1,:,ilat,ilon)*1e3; %g/m2
n=1;
iwp_summit_3_by_3 = dat_ALL.dat_ens(1,:,ilat-n:ilat+n,ilon-n:ilon+n)*1e3; %g/m2
iwp_summit_3_by_3_av = meanNoNan(iwp_summit_3_by_3(:,:,:),3);

model_str = 'u-bz785, wind nudging, abv. BL'; load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_Nudged_u-bz785_all_tot_cloud_amount_in_rad.mat';
dat_ALL=load(load_file,'dat_ens');
totCF_summit = dat_ALL.dat_ens(1,:,ilat,ilon); %
n=1;
totCF_summit_3_by_3 = dat_ALL.dat_ens(1,:,ilat-n:ilat+n,ilon-n:ilon+n); %g/m2
totCF_summit_3_by_3_av = meanNoNan(totCF_summit_3_by_3(:,:,:),3);


%% Save data

file_save = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Greenland_Summit_u-bz785_all_LWP_2-391_timeseries.mat';
save(file_save,'-V7.3','lwp_summit','lwp_summit_3_by_3_av','iwp_summit','iwp_summit_3_by_3_av','totCF_summit','totCF_summit_3_by_3_av','lat_summit','lon_summit','years','months','time_decimal_year');
mat2nc_Dan(file_save,[file_save '.nc']);

%Add units to the file - don't seem to work using ! - copy and run on command line
%ncatted -a units,"lwp_summit_3_by_3_av",c,c,"g/m2 - closest model point to Summit" /home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Greenland_Summit_u-bz785_all_LWP_2-391_timeseries.mat.nc
%ncatted -a units,"lwp_summit",c,c,"g/m2 - average over the 3x3 region centred on Summit" /home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Greenland_Summit_u-bz785_all_LWP_2-391_timeseries.mat.nc

%ncatted -a units,"iwp_summit_3_by_3_av",c,c,"g/m2 - closest model point to Summit" /home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Greenland_Summit_u-bz785_all_LWP_2-391_timeseries.mat.nc
%ncatted -a units,"iwp_summit",c,c,"g/m2 - average over the 3x3 region centred on Summit" /home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Greenland_Summit_u-bz785_all_LWP_2-391_timeseries.mat.nc

%ncatted -a units,"totCF_summit_3_by_3_av",c,c,"fraction between 0 and 1 - closest model point to Summit" /home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Greenland_Summit_u-bz785_all_LWP_2-391_timeseries.mat.nc
%ncatted -a units,"totCF_summit",c,c,"fraction between 0 and 1 - average over the 3x3 region centred on Summit" /home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Greenland_Summit_u-bz785_all_LWP_2-391_timeseries.mat.nc


%% Plot data

figure
set(gcf,'color','w');
plot(time,lwp_summit,'b-');
ylabel('LWP (g m^{-2})');
xlabel('Year');

figure
set(gcf,'color','w');
lwp = reshape(lwp_summit,[12 length(lwp_summit)/12]);
monthly_mean = meanNoNan(lwp,2);
plot(monthly_mean,'bo-');
ylabel('LWP (g m^{-2})');
xlabel('Month');

%%

figure
set(gcf,'color','w');
plot(time,iwp_summit,'b-');
ylabel('IWP (g m^{-2})');
xlabel('Year');

figure
set(gcf,'color','w');
iwp = reshape(iwp_summit,[12 length(iwp_summit)/12]);
monthly_mean = meanNoNan(iwp,2);
plot(monthly_mean,'bo-');
ylabel('IWP (g m^{-2})');
xlabel('Month');

%%

figure
set(gcf,'color','w');
plot(time,totCF_summit,'b-');
ylabel('Total cloud fraction');
xlabel('Year');

figure
set(gcf,'color','w');
totCF = reshape(totCF_summit,[12 length(totCF_summit)/12]);
monthly_mean = meanNoNan(totCF,2);
plot(monthly_mean,'bo-');
ylabel('Total cloud fraction');
xlabel('Month');


