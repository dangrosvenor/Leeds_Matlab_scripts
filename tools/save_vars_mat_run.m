%can use this to save various variables (e.g. following DRIVER_plot_global)

vars_class = 'seasonal data';
vars_class = 'plot global maps data';

tag = 'AMSREvsCF';
tag = 'MODIS_MOD35_15_percent_vsCF';
tag = 'MODIS_MOD06_vsCF';
%tag = 'MODIS_MOD35_vsCF';
tag = 'Nd_multi_year_seasonal_SouthernOcean_L3';
tag = 'Nd_2007_wavelengths';
tag = '_AMSRE_clear_sky';
%tag = 'MODIS_LWP_for_AMSRE_clear_sky';

%tag = 'lowSZA';
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt']; %re
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Nd']; %Nd
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau']; %tau
filename = ['saved_vars_for_Nd_multi_year_seasonal_SouthernOcean_L3'];  %for the 'MODIS Nd multiple year seasonal cycles'
filename = ['saved_vars_for_Nd_2007_all_wavelengths_CF80_seasonal_SouthernOcean_mockL3'];  %
filename = ['saved_vars_for_Nd_2007_all_wavelengths_CF99_CTT_268_seasonal_SouthernOcean_mockL3'];  %
filename = ['saved_vars_Nd2007_all_wlengths_CF80_seasonal_60-70S_60-160W_L3_individ'];  %
filename = ['AMSRE_clear_sky_bias_7percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to
   %considered clear sky.
%filename = ['AMSRE_clear_sky_bias_MODIS_LWP_10percent_VOCALS_region'];  %
filename = ['AMSRE_clear_sky_bias_9percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to
filename = ['AMSRE_clear_sky_bias_5percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to
filename = ['AMSRE_clear_sky_bias_3percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to
filename = ['AMSRE_clear_sky_bias_80percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to
filename = ['AMSRE_clear_sky_bias_extrapolated_VOCALS_region'];  %Extrapolated to zero CF using 5, 7 and 9 % CFs - see amsre_clear_sky_bias_extrapolate.m
   %case of DRIVER_plot_global_maps

filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
filename2 = [filename '_' datestr(now,30) '.mat'];
filename_savevars = [filedir_savevars filename2];

save_or_load = 'save';


ivar=1; clear vars_to_save

% vars_to_save{ivar}='xdat'; ivar=ivar+1;
% vars_to_save{ivar}='ydat'; ivar=ivar+1;
% vars_to_save{ivar}='errordatU'; ivar=ivar+1;
% vars_to_save{ivar}='errordatL'; ivar=ivar+1;
% vars_to_save{ivar}='thresh_str'; ivar=ivar+1;

switch vars_class


    case 'seasonal data'

        vars_to_save{ivar}='season_mean_driver'; ivar=ivar+1;
        vars_to_save{ivar}='season_Ndatap_driver'; ivar=ivar+1;
        vars_to_save{ivar}='season_std_driver'; ivar=ivar+1;
        vars_to_save{ivar}='mean_sig2_driver'; ivar=ivar+1;
        vars_to_save{ivar}='Nspatial_driver'; ivar=ivar+1;
        vars_to_save{ivar}='season_timestd_driver'; ivar=ivar+1;

        vars_to_save{ivar}='season_datatype_driver'; ivar=ivar+1;
        vars_to_save{ivar}='season_vals_timelabs_driver'; ivar=ivar+1;

        vars_to_save{ivar}='set_years'; ivar=ivar+1;
        vars_to_save{ivar}='screening_eval_str'; ivar=ivar+1;
        vars_to_save{ivar}='screen_edits'; ivar=ivar+1;

        vars_to_save{ivar}='thresh_str_driver'; ivar=ivar+1;
        vars_to_save{ivar}='screen_type_driver'; ivar=ivar+1;
        vars_to_save{ivar}='years_required_for_mean_driver'; ivar=ivar+1;

    case 'plot global maps data';
        vars_to_save{ivar}='P_save'; ivar=ivar+1; %save P_save rather than P because 
        %P might get scaled to make new colorbars (inew_cticks)
        vars_to_save{ivar}='Plon2D_edges'; ivar=ivar+1;
        vars_to_save{ivar}='Plat2D_edges'; ivar=ivar+1;
        vars_to_save{ivar}='Plon2D'; ivar=ivar+1;
        vars_to_save{ivar}='Plat2D'; ivar=ivar+1;
        vars_to_save{ivar}='MODIS_varname2_plot'; ivar=ivar+1;
        vars_to_save{ivar}='units_str_plot'; ivar=ivar+1;
        
        vars_to_save{ivar}='cont_dat'; ivar=ivar+1;
        vars_to_save{ivar}='cont_ints'; ivar=ivar+1;
        vars_to_save{ivar}='icontour'; ivar=ivar+1;
        
        vars_to_save{ivar}='Npoints'; ivar=ivar+1;        
        vars_to_save{ivar}='Ndays2'; ivar=ivar+1;                

end

%run the save script.
save_vars_mat