function [time_inds,time_out] = UM_get_time_inds_from_dat_global(dat_global,time_choice,opts)

if exist('opts')
    %Convert all of the variable names in the input structure to actual names
    %for ease of use
    name_struc='opts'; %The name of the structure
    names = eval(['fieldnames(' name_struc ');']);
    for i=1:length(names)
        eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
        eval(eval_str);
    end
    
end


array_in=[]; %just test to get the indices for now.
dim=NaN; %don't need dim if just getting the indices
[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

ierror = 0;
if exist('time_out_check')
    if length(time_out) ~= length(time_out_check)
        ierror=1;
        error_message = ['length(time_out)=' num2str(length(time_out)) ' not equal to length(time_out_check)=' num2str(length(time_out_check))];
    else
        max_diff = maxALL(abs(time_out - time_out_check));
        if ~exist('time_tol') %from opts
            time_tol = 1e-30;
        end
        if max_diff > time_tol
            ierror=1;
            error_message = ['max time diff, max_diff (' num2str(max_diff) ') > time_tol (' num2str(time_tol) ')'];
        end
    end
    if ierror==1
        %message = ['\n***DPG - times do not match the comparison times!***\n'];
        if exist('no_error') & no_error==1
            fprintf(1,error_message);
        else
            error(error_message);
        end
    end
    
end



