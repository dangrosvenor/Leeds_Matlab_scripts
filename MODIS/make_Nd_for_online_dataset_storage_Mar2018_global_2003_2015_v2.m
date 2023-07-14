%This takes the .mat files already produced (using
%make_Nd_for_McCoy_May2017_global_2003_2015.m) and makes them a bit more
%user friendly for online storage.

file_dir_IN = '/home/disk/eos15/d.grosvenor/mock_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/';
file_dir_OUT='/home/disk/eos15/d.grosvenor/mock_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/for_online_storage/';

files = dir([file_dir_IN 'Nd37*.mat']);

for ifile=1:length(files)     
    
    load([file_dir_IN files(ifile).name]);

%[date_str,date_num] = date_from_day_of_year_func(daynum_timeseries3,modisyear_timeseries3);
% Do the time as days since 1st Jan 1970.
%time_vec = datevec(date_num); time_vec=time_vec(:,1:3);
time_vec = time;
date_num = datenum(time_vec);
time = date_num - datenum('01-Jan-1970');
%lat=LAT;
%lon=LON;

%Notes='Produced using make_Nd_for_McCoy_May2017_global_2003_2015.m';
%Notes = '';

%2.1um
%Nd = permute(N_time3,[2 1 3]);
%Nd = permute(Droplet_Number_Concentration.timeseries3,[2 1 3]);
%file_save = [file_dir 'Nd21_' num2str(time(1,1)) '_SZA_LT_65_CF_GT_80_CTH_LT_3.2km.mat'];

%save(file_save,'-V7.3','Nd','time','lon','lat','Nd_1km','Notes','make_mockL3_variables_script_name','multiL2L3_case','multiL2L3_project','daily_averaged_files_loc2');
%mat2nc_Dan(file_save,[file_save '.nc']);

%3.7um
%Nd = permute(N_time3_37,[2 1 3]);
%Nd = permute(Droplet_Number_Concentration_37.timeseries3,[2 1 3]);
Nd = Nd_1km;
file_save = [file_dir_OUT 'Nd37_' num2str(time_vec(1,1)) '_SZA_LT_65_CF_GT_80_CTH_LT_3.2km.mat'];

%save(file_save,'-V7.3','Nd','time','lon','lat','Nd_1km','Notes','make_mockL3_variables_script_name','multiL2L3_case','multiL2L3_project','daily_averaged_files_loc2');
save(file_save,'-V7.3','Nd','time','time_vec','lon','lat');
file_nc = [remove_character(file_save,'.mat','') '.nc'];
mat2nc_Dan(file_save,file_nc);

eval(['!ncatted -a units,''time'',c,c,''days since 1st Jan 1970'' ' file_nc]);
eval(['!ncatted -a units,''Nd'',c,c,''cm^{-3}'' ' file_nc]);
eval(['!ncatted -a units,''time_vec'',c,c,''year month day'' ' file_nc]);

end


