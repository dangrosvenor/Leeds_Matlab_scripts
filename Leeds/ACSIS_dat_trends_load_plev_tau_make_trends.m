function []=ACSIS_dat_trends_load_plev_tau_make_trends()

%% Choose region and land/ocean screening
box_region = '4';
box_region = '18'; %Northern NA as used in ACP paper
%box_region = '19'; %Southern NA as used in ACP paper

land_ocean = 'land+ocean';
land_ocean = 'ocean only';
%land_ocean = 'land only'; 

ACSIS_Robson_paper_choose_regional_box2 %run script - also chooses ylims, etc.

%%
region_str=box_region_str;
var_UM='clisccp';

%Directory to process :-
opts.named_dir =['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/ens_mean/' var_UM '/'];
%Directory to save to
savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_' region_str '_trend_ukesm'];

season='Annual';
dat_PI='';
start_year_dat = 1850;
end_year_dat = 2014;

%% Trend choices

yr_start_trend_box = [1985]; yr_end_trend_box = [2014]; %can select for multiple trends

p_conf = 95; % Confidence limit (%) for the trend significance
nthresh_days = 3;

%%
load_type = 'mat';
load_type = 'merged netCDF';
load_type = 'individual netCDF files'; opts.cat_dim=1;
load_type = 'named netCDF'; %opts.named_file='foo';

savefile = [savefile_pre_str '_' var_UM '.mat'];

clear trends_tau_plev
for itau=0:6
    itau
    for iplev=0:6
        
        opts.named_file=['clisccp_ens_mean_merged_times_p' num2str(iplev) '_tau' num2str(itau) '.nc.nc3'];
        
        nc=netcdf(opts.named_file);        
   
        pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %               
        
        opts.isort_time=1;
        opts.lat_var = 'lat';
        opts.lon_var = 'lon';
        
        opts.time_var='time';
        opts.time_ref = datenum('01-Jan-1850');
        opts.time_fconv = 1; %conversion multiplier to get to days
        dat_global = UM_load_merged_netCDF(opts.named_dir,var_UM,run_type,load_type,[],[],var_UM,opts); %data is ordered [time lat lon]. 180 times (monthly over 15 years)                    
        
        dat_ukesm = dat_global;
        dat_ukesm.years_ukesm_1d = [start_year_dat:end_year_dat];
        dat_ukesm.dat_annual = squeeze(dat_global.dat);        
        gcm_area_UM = dat_global.gcm_area_UM;
           
        ibox=1;
        
        %Run script to choose the data for the selected season and then to average
        %over the box region.
        yr_start_trend_box2 = yr_start_trend_box; yr_end_trend_box2 = yr_end_trend_box;
        ACSIS_Robson_paper_box_means
        
        trends_tau_plev{itau+1,iplev+1} = trend_dat_box{1};
        dat_box_tau_plev{itau+1,iplev+1} = dat_annual_box_ukesm;
            
            
    end

end

save(savefile,'-V7.3'); %save all the variables


