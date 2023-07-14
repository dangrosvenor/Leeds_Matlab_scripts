%Load the data using MODIS_multi_DAY_processL3L2 with make_mockL3_variables_script_name = 'make_mockL3_variables_just_for_Nd_May2017';
% and multiL2L3_case = 'load L3 and concatenate'; for one year.
% Did this manually for each year.

%Already have the Droplet_Number variables (made using 1km tau and reff).
%And 1x1 deg values (made in filtering_data_get) - N_time3 and N_time3_37
%McCoy would prefer Nd21 and 37 in separate files formatted as
%[lon,lat,time] with a time vector of size [N 3] (Y M D) and lat and lon
%(vectors I guess?)

file_dir='/home/disk/eos15/d.grosvenor/mock_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/'

[date_str,date_num] = date_from_day_of_year_func(daynum_timeseries3,modisyear_timeseries3);
time = datevec(date_num); time=time(:,1:3);
lat=LAT;
lon=LON;

Notes='Produced using make_Nd_for_McCoy_May2017_global_2003_2015.m';

%2.1um
Nd = permute(N_time3,[2 1 3]);
Nd_1km = permute(Droplet_Number_Concentration.timeseries3,[2 1 3]);
file_save = [file_dir 'Nd21_' num2str(time(1,1)) '_SZA_LT_65_CF_GT_80_CTH_LT_3.2km.mat'];

save(file_save,'-V7.3','Nd','time','lon','lat','Nd_1km','Notes','make_mockL3_variables_script_name','multiL2L3_case','multiL2L3_project','daily_averaged_files_loc2');
mat2nc_Dan(file_save,[file_save '.nc']);

%3.7um
Nd = permute(N_time3_37,[2 1 3]);
Nd_1km = permute(Droplet_Number_Concentration_37.timeseries3,[2 1 3]);
file_save = [file_dir 'Nd37_' num2str(time(1,1)) '_SZA_LT_65_CF_GT_80_CTH_LT_3.2km.mat'];

save(file_save,'-V7.3','Nd','time','lon','lat','Nd_1km','Notes','make_mockL3_variables_script_name','multiL2L3_case','multiL2L3_project','daily_averaged_files_loc2');
mat2nc_Dan(file_save,[file_save '.nc']);


