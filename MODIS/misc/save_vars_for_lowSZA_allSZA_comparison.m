%save some variables for a comparison between <65 deg SZA data (new L3) and data at
%all SZA (as for normal L3)

save_dir='/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/';
save_file=[save_dir 'low_sza_vs_all_sza_comparison_terra_2007_SZA_Paper_screenings_LOW.mat']; data_tag = '_lowSZA';
%save_file=[save_dir 'low_sza_vs_all_sza_comparison_terra_2007_SZA_Paper_screenings_ALL.mat']; data_tag = '_allSZA';



save_vars_for_lowSZA_allSZA_comparison_varlist

Nd_16 = N_time3_16;
Nd_21 = N_time3;
Nd_37 = N_time3_37;

%make tau = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
%slightly different NaN values)
Nd_nan = find(isnan(Nd_16)==1);
Nd_21(Nd_nan) = NaN;
Nd_37(Nd_nan) = NaN;
Nd_nan = find(isnan(Nd_21)==1);
Nd_16(Nd_nan) = NaN;
Nd_37(Nd_nan) = NaN;
Nd_nan = find(isnan(Nd_37)==1);
Nd_16(Nd_nan) = NaN;
Nd_21(Nd_nan) = NaN;


Re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3;
Re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3;
Re_37 = Cloud_Effective_Radius_37_Liquid_Mean.timeseries3;

%make tau = NaN when we have NaNs in re_1.6 and re_3.7 (re_2.1 has
%slightly different NaN values)
re_nan = find(isnan(Re_16)==1);
Re_21(re_nan) = NaN;
Re_37(re_nan) = NaN;
re_nan = find(isnan(Re_21)==1);
Re_16(re_nan) = NaN;
Re_37(re_nan) = NaN;
re_nan = find(isnan(Re_37)==1);
Re_16(re_nan) = NaN;
Re_21(re_nan) = NaN;



app_str=')';
for ivar=1:length(save_vars)
%    eval(['save(''' save_file ''',''' save_vars{ivar}  ''',''-V7.3''' app_str]);   
    eval(['save_vars_for_lowSZA_allSZA_comparison_func(save_file,' save_vars{ivar} ',''' save_vars{ivar} ''',data_tag,app_str);']);
    app_str=',''-APPEND'')';
end
    
disp('*** Done save ***');



