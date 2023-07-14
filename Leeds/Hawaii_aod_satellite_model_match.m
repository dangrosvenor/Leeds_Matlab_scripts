%% AOD MODIS matching
% Load the modis AOD data from the .mat files
% Load the MODIS AOD data
filedir='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/Hawaii_Dec2020_ADVANCE/';
%file_name_h5='MYD04_L2.A2020356.0055.061.2020357191130.hdf';
save_name = 'Hawaii_all_L2_04_data.mat';
%modis_aod_file = [filedir save_name];
modis_aod = load([filedir save_name]);

% Load the model data from the .mat files for PI
UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';
%um_case_PI = 'u-cj086';  %volcano OFF   
um_case_PI = 'u-ch765';  %volcano OFF, orog ON
% Load the model AOD data from the .mat file
aod_save_file = [UM_base_dir um_case_PI '/aod550_total.mat'];
model_aod_PI = load(aod_save_file);

% Load the model data from the .mat files for PD
UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';
%um_case_PD = 'u-co295';  %volcano ON    |
um_case_PD = 'u-cr138';  %volcano ON, orog ON    |
% Load the model AOD data from the .mat file
aod_save_file = [UM_base_dir um_case_PD '/aod550_total.mat'];
model_aod_PD = load(aod_save_file);

var_str='aod';

% Match the PI data for different coarse grain values
coarse_str = '3x3';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_aod_PI,UM_base_dir,um_case_PI,coarse_str);
coarse_str = '2x2';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_aod_PI,UM_base_dir,um_case_PI,coarse_str);

% Match the PD data for different coarse grain values
coarse_str = '3x3';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_aod_PD,UM_base_dir,um_case_PD,coarse_str);
coarse_str = '2x2';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_aod_PD,UM_base_dir,um_case_PD,coarse_str);


%% SO2 MODIS matching
% Load the model data from the .mat files for PI
UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';
%um_case_PI = 'u-cj086';  %volcano OFF  
um_case_PI = 'u-ch765';  %volcano OFF, orog ON
% Load the model AOD data from the .mat file
so2_save_file = [UM_base_dir um_case_PI '/so2.mat'];
model_so2_PI = load(so2_save_file);

% Load the model data from the .mat files for PD
UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';
%um_case_PD = 'u-co295';  %volcano ON   
um_case_PD = 'u-cr138';  %volcano ON, orog ON    
% Load the model AOD data from the .mat file
so2_save_file = [UM_base_dir um_case_PD '/so2.mat'];
model_so2_PD = load(so2_save_file);

var_str='so2';

% Match the PI data for different coarse grain values
coarse_str = '3x3';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_so2_PI,UM_base_dir,um_case_PI,coarse_str);
coarse_str = '2x2';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_so2_PI,UM_base_dir,um_case_PI,coarse_str);

% Match the PD data for different coarse grain values
coarse_str = '3x3';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_so2_PD,UM_base_dir,um_case_PD,coarse_str);
coarse_str = '2x2';
[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_so2_PD,UM_base_dir,um_case_PD,coarse_str);



