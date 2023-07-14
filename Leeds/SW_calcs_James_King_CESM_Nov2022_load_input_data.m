% SW_calcs_James_King_CESM_Nov2022_load_input_data.m
% Run from SW_calcs_James_King_Nov2022_RUN.m

%data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/James_King_CESM/';
clear optionals
optionals.iget_time_only=1;

%model_type = 'Max forest ssp3 2050';
model_type = 'Max forest ssp3 2095';

switch model_type
    case 'Max forest ssp3 2050'
        
        %Run both once to get the time indices
        year='2050';
        %year='2095';
        %UM_run_PI = ['ssp3_base_' year]; file_suffix_PI = '_3-5yrs'; %control run
        UM_run_PI = ['ssp3_base_' year]; file_suffix_PI = ''; %control run
        file_prefix_PI = ['CAM_chem_SSP3-7.0_' year '_ssp_landuse_all_forcings_SP_mode.'];
        [nT01] = SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02(UM_run_PI,file_prefix_PI,file_suffix_PI,optionals);
        CESM_run_PI = UM_run_PI;
        
        %UM_run_PD = ['maxforest_adjusted_ssp3_' year]; file_suffix_PD = '_3-5yrs'; %control run
        UM_run_PD = ['maxforest_adjusted_ssp3_' year]; file_suffix_PD = ''; %control run
        file_prefix_PD = ['CAM_chem_SSP3-7.0_' year '_maxforest_complete_adjust_SP_mode.'];
        [nT02] = SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02(UM_run_PD,file_prefix_PD,file_suffix_PD,optionals);
        CESM_run_PD = UM_run_PD;
        
    case 'Max forest ssp3 2095'
        
        %Run both once to get the time indices        
        year='2095';
        %UM_run_PI = ['ssp3_base_' year]; file_suffix_PI = '_3-5yrs'; %control run
        UM_run_PI = ['ssp3_base_' year]; file_suffix_PI = ''; %control run
        file_prefix_PI = ['CAM_chem_SSP3-7.0_' year '_ssp_landuse_all_forcings_SP_mode.'];
        [nT01] = SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02(UM_run_PI,file_prefix_PI,file_suffix_PI,optionals);
        CESM_run_PI = UM_run_PI;
        
        %UM_run_PD = ['maxforest_adjusted_ssp3_' year]; file_suffix_PD = '_3-5yrs'; %control run
        UM_run_PD = ['maxforest_adjusted_ssp3_' year]; file_suffix_PD = ''; %control run
        file_prefix_PD = ['CAM_chem_SSP3-7.0_' year '_maxforest_complete_adjust_SP_mode.'];
        [nT02] = SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02(UM_run_PD,file_prefix_PD,file_suffix_PD,optionals);
        CESM_run_PD = UM_run_PD;        
        
        
end

%Run again with the time indices limited.
optionals.iget_time_only = 0;
optionals.nT = min(nT01.nT,nT02.nT);

[dat_PI] = SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02(UM_run_PI,file_prefix_PI,file_suffix_PI,optionals);
[dat_PD] = SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02(UM_run_PD,file_prefix_PD,file_suffix_PD,optionals);

 
i_Nd_plot = 0;
if i_Nd_plot == 1
   lats = [18 60];
   lons = [-75 0];
   var = 'Nd';   
   [me_CESM] = 1e-6*mean_area_weighted_timeseries_model_data(dat_PI,var,lats,lons); me_CESM2 = meanNoNan(me_CESM,2);
   [me_UM] = mean_area_weighted_timeseries_model_data(dat_UM_PI,var,lats,lons); me_UM2 = meanNoNan(me_UM,2);
   
   i=1; clear leg_str
   figure('color','w')
   increase_font_size_map_figures
   plot(me_CESM,'linewidth',3); leg_str{i}=['CESM, \mu=' num2str(me_CESM2)]; i=i+1;   
   hold on
   plot(me_UM,'r','linewidth',3); leg_str{i}=['UKESM, \mu=' num2str(me_UM2)]; i=i+1
   title([num2str(lats(1)) ' to ' num2str(lats(2)) 'N, ' num2str(lons(1)) ' to ' num2str(lons(2)) 'E']);
   xlabel('Time');
   ylabel('N_d (cm^{-3})');
   legend(leg_str);
    
end
