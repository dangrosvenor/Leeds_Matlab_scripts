function [model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func01(var_str,modis_aod,model_aod,UM_base_dir,um_case_PD,coarse_str)

%% Load the model and modis data from the .mat files
%UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';
%um_case_PD = 'u-co295';  %volcano ON    

% Load the model AOD data from the .mat file
%aod_save_file = [UM_base_dir um_case_PD '/aod550_total.mat'];
%model_aod = load(aod_save_file);

% Load the MODIS AOD data
%filedir='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/Hawaii_Dec2020_ADVANCE/'    
%file_name_h5='MYD04_L2.A2020356.0055.061.2020357191130.hdf';
%save_name = 'Hawaii_all_L2_04_data.mat';
%modis_aod = load([filedir save_name]);

eval_str=['model_lat2D = model_aod.gcm_Plat2D_' coarse_str ';']; eval(eval_str);
eval_str=['model_lon2D = model_aod.gcm_Plon2D_' coarse_str ';']; eval(eval_str);
eval_str=['model_dat = model_aod.dat_' coarse_str ';']; eval(eval_str);


[model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func02(modis_aod,model_lat2D,model_lon2D,model_dat,model_aod.time_out_var);

model_save_file = [UM_base_dir um_case_PD '/modis_' var_str '_matches_' coarse_str '.mat'];
fprintf(1,['\nSaving to file: ' model_save_file '\n']);
save(model_save_file,'-V7.3','model_match','modis_match','d_match');


