%this is usually run using multi_run_gcm_time_chunk



gcm_saved_filename_file = [filedir_gcm_load 'y' file_year file_month '_timechunks_' am3_dataset '_' datestr(now,30) '.mat'];
%may want to save this filename if processing multiple half-year files

gcm_time_read = nc_inst{'time'}(:); %

gcm_str = chunk_model;

%if nT==1
%    disp('Don''t choose nT==1 - read_am3T_savemem doens''t work in that case (or fix this!)')
%    return
%end

nTall = length(gcm_time_read);

nstarts = [1:nT:nTall]; %nstarts(end) will not necessarily be =nTall
nstarts(end+1)=nTall+1; %since each time index vector will run to t(i+1)-1




%% the time loop for loading, processing and saving the time chunks one at
%% a time
if ~exist('ioverride_multi_chunk') | ioverride_multi_chunk==0
    itam32=1; %otherwise keep incrementing for the whole year
end

for itam3=1:length(nstarts)-1

    
    fprintf(1,'\nProcessing chunk %d of %d....',itam3,length(nstarts)-1);
    
    gcm_idays=[nstarts(itam3):nstarts(itam3+1)-1];

    %read in the data for these indices
    ioverride_read_am3=1;
    
    if length(findstr(chunk_model,'AM3'))>0
        if cosp_flag==1
            read_AM3T_cosp_savemem
        else
            iread_cosp=0;
            read_AM3T_cosp_savemem
        end
    elseif length(findstr(chunk_model,'CAM'))>0
        if cosp_flag==1
            read_cosp_cam
        else
            read_cam
        end
    else
        error('*** Need to specify the reading function for this model in am3_multi_load_time_chunks ***');
    end
    
%     switch chunk_model
%         case {'AM3','AM3_CLUBBv2_COSP_100km','AM3_100km_20130110'} %these had COSP turned on
%             %iread_cosp=1;
%             read_AM3T_cosp_savemem
%         case {'AM3_CLUBB','AM3_CLUBBv2','AM3_CLUBBv1_1deg','AM3_CLUBBv1_2deg'} %the pre-July 2013 CLUBB runs had no COSP
%             iread_cosp=0;
%             read_AM3T_cosp_savemem       
%         case {'CAM5_COSP','CAM5_CLUBB_COSP','CAM5_CLUBBv2_COSP'}
%             read_cosp_cam
%         case {'CAM5','CAM5_CLUBB'}
%             read_cam
%         otherwise
%             error('*** Need to specify the reading function for this model in am3_multi_load_time_chunks ***');            
%     end

    
    %process the data
    gcm_process
    
    %save the data
    ioverride_saveload=1;
    if exist('h_half')
            save_gcm_process_vars_flag=1; %=1 means that it saves h_half, isscp_lwp, cloud thicknesses etc    
    else
            save_gcm_process_vars_flag=2; %=2 means that it doesn't save h_half, isscp_lwp, cloud thicknesses etc 
            %since they aren't available.
    end
    load_gcm_process_vars_flag=0; %also need to make sure we set the load flag too.
    save_gcm_processed_data
    
    %store the filenames
    
    gcm_year_cases{itam32} = ['y' file_year]; gcm_month_cases{itam32} = file_month; gcm_model_cases{itam32} = chunk_model;
    itam32=itam32+1; %
    chunk_name = [gcm_year_cases{itam32-1} '_' gcm_model_cases{itam32-1} '_' gcm_month_cases{itam32-1} '_chunk_' num2str(itam32-1)];    
    
    %save the filepath to the .mat file
    eval(['am3_filepaths_chunk.' chunk_name ' = savename_gcm']);
end





