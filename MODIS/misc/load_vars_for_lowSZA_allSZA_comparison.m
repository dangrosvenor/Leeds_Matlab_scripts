%save_vars_for_lowSZA_allSZA_comparison.m used to save the variables.


%list of the variables
save_vars_for_lowSZA_allSZA_comparison_varlist

data_tags = {'_lowSZA','_allSZA'};
save_dir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/';
%with just CF>0.8 screening
save_files={'low_sza_vs_all_sza_comparison_terra_2007_LOW.mat','low_sza_vs_all_sza_comparison_terra_2007_ALL.mat'};
%with SZA_paper screenings
save_files={'low_sza_vs_all_sza_comparison_terra_2007_SZA_Paper_screenings_LOW.mat','low_sza_vs_all_sza_comparison_terra_2007_SZA_Paper_screenings_ALL.mat'};



for itag=1:length(data_tags)

    save_file = [save_dir save_files{itag}];
    app_str=')';
    for ivar=1:length(save_vars)
            eval(['load(''' save_file ''',''' save_vars{ivar}  data_tags{itag} ''');']);
    end

end
    

eval(['MLAT=MLAT' data_tags{1} ';']);
eval(['MLON=MLON' data_tags{1} ';']);
eval(['daynum_timeseries3=daynum_timeseries3' data_tags{1} ';']);
eval(['modisyear_timeseries3=modisyear_timeseries3' data_tags{1} ';']);
daynum_timeseries3_MODIS = daynum_timeseries3;
LAT=MLAT;
LON=MLON;

disp('*** Done load ***');