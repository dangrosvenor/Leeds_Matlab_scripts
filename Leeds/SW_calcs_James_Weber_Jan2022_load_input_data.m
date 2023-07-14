% SW_calcs_James_Weber_Jan2022_load_input_data.m
% Run from SW_calcs_James_Weber_Jan2022_RUN.m





model_type = 'Strat-Trop';
model_type = 'CRI-Strat 2';
model_type = 'Max forest ssp3 2050';
%model_type = 'Max forest ssp3 2095';
model_type = 'Max forest ssp126 2050';
model_type = 'Max forest ssp126 2095';

switch model_type
    case 'Strat-Trop'        
        data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/';
        
        UM_run = 'u-ck596'; file_suffix='_8_years'; %control run
        dat_PI = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PI = UM_run;
        
        UM_run = 'u-ck597'; file_suffix='_8_years'; %perturbed run
        dat_PD = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PD = UM_run;
        
    case 'CRI-Strat 2'
        data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/';
        
        % CRI-Strat 2 runs
        %UM_run = 'u-ce358'; file_suffix=''; %control run
        UM_run = 'u-ck598'; file_suffix='_8_years'; %control run
        dat_PI = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PI = UM_run;
        
        %UM_run = 'u-ce419'; file_suffix=''; %perturbed run
        UM_run = 'u-ck599'; file_suffix='_8_years'; %perturbed run
        dat_PD = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PD = UM_run;
        
    case 'Max forest ssp3 2050'
        data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/James_King_CESM/forests_UKESM/';
        
        % CRI-Strat 2 runs
        %UM_run = 'u-ce358'; file_suffix=''; %control run
        UM_run = 'u-co022'; file_suffix=['_2052_2064_' UM_run(3:end)]; %control run
        dat_PI = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PI = UM_run;
        
        %UM_run = 'u-ce419'; file_suffix=''; %perturbed run
        UM_run = 'u-cr012'; file_suffix=['_2052_2064_' UM_run(3:end)]; %perturbed run
        dat_PD = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PD = UM_run;  
        
    case 'Max forest ssp3 2095'
        data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/James_King_CESM/forests_UKESM/';
                        
        UM_run = 'u-co023'; file_suffix=['_2097_2109_' UM_run(3:end)]; %control run
        dat_PI = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PI = UM_run;
                
        UM_run = 'u-cr013'; file_suffix=['_2097_2109_' UM_run(3:end)]; %perturbed run
        dat_PD = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PD = UM_run; 
        
 case 'Max forest ssp126 2050'
        data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/James_King_CESM/forests_UKESM/';
        
        % CRI-Strat 2 runs
        %UM_run = 'u-ce358'; file_suffix=''; %control run
        UM_run = 'u-co112'; file_suffix=['_2052_2064_' UM_run(3:end)]; %control run
        dat_PI = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PI = UM_run;
        
        %UM_run = 'u-ce419'; file_suffix=''; %perturbed run
        UM_run = 'u-cr140'; file_suffix=['_2052_2064_' UM_run(3:end)]; %perturbed run
        dat_PD = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PD = UM_run;  
        
    case 'Max forest ssp126 2095'
        data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/James_King_CESM/forests_UKESM/';
                        
        UM_run = 'u-co113'; file_suffix=['_2097_2109_' UM_run(3:end)]; %control run
        dat_PI = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PI = UM_run;
                
        UM_run = 'u-cr141'; file_suffix=['_2097_2109_' UM_run(3:end)]; %perturbed run
        dat_PD = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir);
        UM_run_PD = UM_run;         
end