%save processed GCM data from multiple years
%variables to save are listed in gcm_list_of_saved_variables

try

if ~exist('ioverride_saveload') | ioverride_saveload==0
    load_gcm_process_vars_flag=0; %whether have ran gcm_process or not (as whether to save those vars or not)
    save_gcm_process_vars_flag=2;
    cosp_flag=0;
end


%selects files for load_saved_modis_vars.m & modis_make_monthly_averages_multi_year.m
% 

savedir_gcm = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/';
savename_gcm = ['saved_gcm_processed_data_' gcm_str '_' am3_dataset '_' nc_inst_file '_' datestr(now,30) '.mat']
savefilename_gcm = [savedir_gcm savename_gcm];

% -------------------------------------------------------------------------
% script that lists the names of all the variables to save in
% gcm_var{i}
gcm_list_of_saved_variables
% -------------------------------------------------------------------------

for ivar=1:length(gcm_var)
    fprintf(1,['\nSaving variable ' num2str(ivar) ' of ' num2str(length(gcm_var))]);


    if ivar==1
        eval_str = ['save(''' savefilename_gcm ''',''' gcm_var{ivar} ''',''-V7.3'');'];
    else
        eval_str = ['save(''' savefilename_gcm ''',''' gcm_var{ivar} ''',''-APPEND'',''-V7.3'');'];
    end
    
    eval(eval_str);


end

fprintf(1,'\nDone\n');

catch save_error
    clear ioverride_saveload
    rethrow(save_error);
end

