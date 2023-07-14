load('name of .mat file for the modelrun in questoin');
%this will contain a gcm_saved_filename_file_YYYY for each year
%Now will just have to add one filename for each model run - can write this
%to screen at the end, but also to a log file

for iyears_load_model_run=1:length(years_load_model_run);

    eval(['load_file{iload_file} = gcm_saved_filename_file_' num2str(years_load_model_run(iyears_load_model_run)) ';'])
    iload_file=iload_file+1;
end