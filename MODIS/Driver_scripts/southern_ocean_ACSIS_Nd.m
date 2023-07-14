thresh_LAT=[-70 -30]; thresh_LON=[-1e9 1e9];


%Calculate the monthly means
[mon_me,years]=calc_monthly_means(Nd_PD_ALL,time_out);
[mon_me_SO] = monthly_means_restrict_lat_lon(years,mon_me,gcm_Plat2D_UM,gcm_Plon2D_UM,thresh_LAT,thresh_LON);

%regrid the model data to the MODIS grid for filtering out later of regions
%with little obs data
clear mon_me_model_MODIS_grid
for iy=1:length(years)  
    for im=1:12
        mon_me_model_MODIS_grid{iy}(:,:,im) = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,mon_me{iy}(:,:,im),gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE);
    end
end
               
figure
dat = [mon_me_SO{2}(1:3) mon_me_SO{1}(4:end)];
plot(dat,'ko-');

figure
dat = [mon_me_SO{1} mon_me_SO{2} mon_me_SO{3} mon_me_SO{4} ]/1e6; %per cc
plot(dat,'ko-');



%% Now deal with the obs data - make into monthly means.
dat_modis = modis_loaded.Droplet_Number_Concentration_37.timeseries3; %CF>80
%dat_modis = modis_loaded_CF0.Droplet_Number_Concentration_37.timeseries3;
[date_str,date_num] = date_from_day_of_year_func(modis_loaded_CF0.daynum_timeseries3,modis_loaded.modisyear_timeseries3);
[mon_me_MODIS,years_MODIS]=calc_monthly_means(dat_modis,date_num);
[mon_me_SO_MODIS] = monthly_means_restrict_lat_lon(years_MODIS,mon_me_MODIS,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,thresh_LAT,thresh_LON);


    

%% Using the data given to Jane - screened for sea-ice etc.
str_2137='21';
%str_2137='37';
nthresh_days = 5;

file_dir='/home/disk/eos1/d.grosvenor/mock_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/';


mon_me_filter = mon_me_model_MODIS_grid;

years_MODIS2=[2009 2010 2011 2012];
clear mon_me_MODIS2 mon_me_MODIS2_Ndatap
for iy=1:length(years_MODIS2)
    year_str = num2str(years_MODIS2(iy));
    filename = [file_dir 'Nd_monthly_' str_2137 '_1deg_' year_str '_SZA_LT_65_CF_GT_80_CTH_LT_3.2km_screened_for_seaice__2week_max.mat.nc'];
    nc=netcdf(filename);
    mon_me_MODIS2{iy} = nc{'Nd_1deg_mean'}(:);
    mon_me_MODIS2_Ndatap{iy} = nc{'Nd_1deg_Ndatap'}(:);
    inan = find(mon_me_MODIS2_Ndatap{iy} < nthresh_days);
    mon_me_MODIS2{iy}(inan)=NaN;
    mon_me_filter{iy}(inan)=NaN;
end

[mon_me_SO_MODIS2] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,thresh_LAT,thresh_LON);
[mon_me_SO_MODIS2_Ndatap] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2_Ndatap,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,thresh_LAT,thresh_LON);
[mon_me_filter_SO] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_filter,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,thresh_LAT,thresh_LON);

%% figure with MODIS and model - WITH model filtering for MODIS locations.
nthresh_days_monthly = 3;
nthresh_days_monthly = 2;
%nthresh_days_monthly = 1;
%nthresh_days_monthly = 0;


figure

dat_model = mon_me_filter_SO;
%dat_model = mon_me_SO;

dat=[];
for i=1:length(dat_model)
    dat = [dat dat_model{i}];
end
dat=dat/1e6;


dat_MODIS=[];
ndat_MODIS=[];
for i=1:length(mon_me_SO_MODIS2)
    dat_MODIS = [dat_MODIS mon_me_SO_MODIS2{i} ];
    ndat_MODIS = [ndat_MODIS mon_me_SO_MODIS2_Ndatap{i} ];
end
dat_MODIS(ndat_MODIS<nthresh_days_monthly) = NaN;
dat(ndat_MODIS<nthresh_days_monthly) = NaN;
plot(dat,'ko-');
hold on
plot(dat_MODIS,'rx-');
grid on

xlabel('Month');
ylabel('Nd (cm^{-3})');
clear leg_str
leg_str{1}='Model';
leg_str{2}='MODIS 3.7um';
title(['Model filtered for same monthly locations as MODIS, nthresh days monthly=' num2str(nthresh_days_monthly) ', ' str_2137 ' MODIS data']);

set(gca,'xlim',[0 50]);

%% figure with MODIS and model - no model filtering for MODIS locations.
nthresh_days_monthly = 2;
%nthresh_days_monthly = 1;
%nthresh_days_monthly = 0;


figure

%dat_model = mon_me_filter_SO;
dat_model = mon_me_SO;

dat=[];
for i=1:length(dat_model)
    dat = [dat dat_model{i}];
end
dat=dat/1e6;


dat_MODIS=[];
ndat_MODIS=[];
for i=1:length(mon_me_SO_MODIS2)
    dat_MODIS = [dat_MODIS mon_me_SO_MODIS2{i} ];
    ndat_MODIS = [ndat_MODIS mon_me_SO_MODIS2_Ndatap{i} ];
end
dat_MODIS(ndat_MODIS<nthresh_days_monthly) = NaN;
dat(ndat_MODIS<nthresh_days_monthly) = NaN;
plot(dat,'ko-');
hold on
plot(dat_MODIS,'rx-');
grid on

xlabel('Month');
ylabel('Nd (cm^{-3})');
clear leg_str
leg_str{1}='Model';
leg_str{2}='MODIS 3.7um';
title(['Model NOT filtered for same monthly locations as MODIS, nthresh days monthly=' num2str(nthresh_days_monthly) ', ' str_2137 ' MODIS data']);

set(gca,'xlim',[0 50]);


