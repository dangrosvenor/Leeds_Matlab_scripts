%Based on UM_Iceland_plot_various_based_on_ACSIS_script.m

%Useful parts of the code
% - time mean variables and deltas and 1d PDFs of those - search for :-
%   "Time mean dLWP, also filtering for LWP"

%Pre-process using
%   UM_quick_plot_global.m

%Add variables to climits are set here too :-
%   UM_var_defs.m

%runs this plotting script - lat/lon for map also set here. FOR PDFs, etc. it is set below :-
  %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
  
% Here, PI will refer to volcano OFF, PD volcano ON

UM_ACSIS_SW_vs_cloud_properties_global_DEFAULTS

%savedir_date=['/home/disk/eos15/d.grosvenor/UM/Hawaii/plots_' datestr(now,30) '/'];
%eval(['!mkdir ' savedir_date]);

savedir_date = '/home/disk/eos15/d.grosvenor/UM/Hawaii/plots_20211213T070655/';
UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';

i_Iceland=1;
icoarse_grain=0;
M_coarse_grain=4; N_coarse_grain=4;
iplot_mgrid_lines_DRIVER=1;
ioverride_ticks_DRIVER=0;
%proj_type_DRIVER; 
  
i_select_region=0; %whether to pick out a sub-domain from the model domain - however, seems to run out of memory doing this...
icloud_states=0; %whether to run the code looking at clodu states, etc.
isave_Dan_Jones=0;
SW_surf_or_TOA = 'TOA'; %choose whether to look at the surfcae SW or the TOA

load_type = 'mat';
load_type = 'merged netCDF';

cloud_input = 'UM'; %Use the usual UM low, mid and high cloud fractions.
%cloud_input = 'CALIPSO'; %Use the UM COSP CALIPSO values
%cloud_input = 'MODIS'; %Use the COSP MODIS values

%[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

amsre_data = 0;

ical_data=0; %Whether are using CALIPSO data - load using script :- read_calipso_monthly_night_IPSL.m for average values
    %first. This puts the data into these fields :- 
    % cllcalipso_monthly_AVERAGE, clmcalipso_monthly_AVERAGE, clhcalipso_monthly_AVERAGE
    % Is average of day and night values for now (have separate ones too).

iceres_data=0; %whether to load CERES data

% Here, PI will refer to volcano OFF, PD volcano ON

%Choose a variable to use to get the times - could be changed below
%var_UM_DRIVER = 'LWP';
%var_UM_DRIVER = 'accum_number_ukca';
%var_UM_DRIVER = 'SW_down_surf';
%var_UM_DRIVER = 'Nd_lwc_weighted_UKCA';
var_UM_DRIVER = 'LWP_time_only'; %use the LWP file, but only for loading the time
var_UM_DRIVER = 'LS_surf_rain_rate_time_only'; %use the file, but only for loading the time


%% Choose runs
run_set = 'global runs Sep 2018';
run_set = 'nested UKCA runs Sep 2018'; %u-ba333 and u-ba050
%run_set = 'nested UKCA-CASIM runs Sep 2018';
%run_set = 'nested UKCA-CASIM-fixed_sed runs Nov 2018'; %u-bc309 and u-bc308
%run_set = 'nested UKCA vs CASIM, volc OFF runs Sep 2018';
%run_set = 'nested UKCA vs CASIM, volc ON runs Sep 2018';
run_set = 'Blending option=3, volcano starting 12 UTC'; %w/orog low emission height
%run_set = 'Blending option=3, volcano starting 12 UTC, no orography'; %low emission height
%run_set = 'Blending option=3, volcano starting 12 UTC, no orography, adjusted emission height'; %higher altitude emissions
%run_set = 'Blending option=3, volcano starting 12 UTC, adjusted emission height'; %Raised emission heights
%run_set = 'Blending option=3, volcano starting 12 UTC, no orography, no direct effect'; %low emission height
run_set = 'Blending option=3, new volcano start time, reduced emissions'; %w/orog low emission height, reduce emission rates
run_set = 'Blending option=3, new volcano start time, reduced emissions further'; %w/orog low emission height, reduce emission rates
run_set = 'Blending option=3, new volcano start time, reduced emissions further, orog OFF'; %NO orog low emission height, reduce emission rates, orog ON.
run_set = '1st June runs'; %As for above, but 1st June runs.


Nd_var_str='UKCA';
Nd_var_str2='';

switch run_set
    case 'global runs Sep 2018'
        um_case_PI = 'u-ba333_glm'; run_type_DRIVER = 'global'; %volcano OFF
        um_case_PD = 'u-ba050_glm'; run_type_DRIVER = 'global'; %volcano ON
    case 'nested UKCA runs Sep 2018'
        um_case_PI = 'u-ba333'; run_type_DRIVER = 'nested'; %volcano OFF
        um_case_PD = 'u-ba050'; run_type_DRIVER = 'nested'; %volcano ON
    case 'nested UKCA-CASIM runs Sep 2018'
        um_case_PI = 'u-ba458'; run_type_DRIVER = 'nested'; %volcano OFF
        um_case_PD = 'u-ba335'; run_type_DRIVER = 'nested'; %volcano ON
    case 'nested UKCA vs CASIM, volc OFF runs Sep 2018'
        um_case_PI = 'u-ba333'; run_type_DRIVER = 'nested'; %UKCA volcano OFF
        um_case_PD = 'u-ba458'; run_type_DRIVER = 'nested'; %volcano ON
    case 'nested UKCA vs CASIM, volc ON runs Sep 2018'
        um_case_PI = 'u-ba050'; run_type_DRIVER = 'nested'; %UKCA volcano OFF
        um_case_PD = 'u-ba335'; run_type_DRIVER = 'nested'; %volcano ON
    case 'nested UKCA-CASIM-fixed_sed runs Nov 2018'
        um_case_PI = 'u-bc309'; run_type_DRIVER = 'nested'; %volcano OFF
        um_case_PD = 'u-bc308'; run_type_DRIVER = 'nested'; %volcano ON   
        Nd_var_str='CASIM';
        
    case 'Blending option=3, volcano starting 12 UTC'
        um_case_PI = 'u-ch765'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF
        um_case_PD = 'u-ch764'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top';
        
        volc_ON_no_orog = 'u-cj085'; %800m emisison height
        volc_OFF_no_orog = 'u-cj086'; 
        %volc_ON_no_orog = 'u-cj725'; %800+1151m
        
        orog_str = 'orography ON'
        emission_str = 'low altitude emissions';
        
    case 'Blending option=3, volcano starting 12 UTC, no orography'
        um_case_PI = 'u-cj086'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF
        um_case_PD = 'u-cj085'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'; 
        
        volc_OFF_no_orog = 'u-cj086';
        volc_ON_no_orog = 'u-cj085';
        
        orog_str = 'orography OFF'
        emission_str = 'low altitude emissions';
        
    case 'Blending option=3, volcano starting 12 UTC, no orography, adjusted emission height'
        um_case_PI = 'u-cj086'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF
        um_case_PD = 'u-cj725'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top';  
        
        orog_str = 'orography OFF'
        emission_str = 'high altitude emissions';
        
    case 'Blending option=3, volcano starting 12 UTC, adjusted emission height'
        um_case_PI = 'u-ch765'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF
        um_case_PD = 'u-cl529'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';
        Nd_var_str2='_total_column_to_zdomain_top';
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top';
        
        %volc_ON_no_orog = 'u-cj085'; %800m emisison height
        volc_OFF_no_orog = 'u-cj086'; 
        volc_ON_no_orog = 'u-cj725'; %800+1151m
        
        orog_str = 'orography ON'
        emission_str = 'high altitude emissions';    
        
    case 'Blending option=3, volcano starting 12 UTC, no orography, no direct effect'
        um_case_PI = 'u-cn100'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF
        um_case_PD = 'u-cn099'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'; 
        
        volc_OFF_no_orog = 'u-cn100';
        volc_ON_no_orog = 'u-cn099';
        
        orog_str = 'orography OFF'
        emission_str = 'low altitude emissions';   
        
    case 'Blending option=3, new volcano start time, reduced emissions'
        um_case_PI = 'u-ch765'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF, orog ON
        um_case_PD = 'u-co295'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON. orog ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'; 
        
        volc_OFF_no_orog = 'u-cj086';
        volc_ON_no_orog = 'u-co296'; %N.B. - not run yet.
        
        orog_str = 'orography OFF'
        emission_str = 'low altitude emissions'; 
        
    case 'Blending option=3, new volcano start time, reduced emissions further'
        um_case_PI = 'u-ch765'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF, orog ON
        um_case_PD = 'u-cr138'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON. orog ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'; 
        
        volc_OFF_no_orog = 'u-cj086';
        volc_ON_no_orog = 'u-cr139'; %N.B. - not run yet so leave these
        %out for now.
        
        
        orog_str = 'orography OFF'
        emission_str = 'low altitude emissions'; 
        
    case 'Blending option=3, new volcano start time, reduced emissions further, orog OFF'
        um_case_PI = 'u-cj086'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF, orog ON
        um_case_PD = 'u-cr139'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON. orog ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'; 
        
        %volc_OFF_no_orog = 'u-cj086';
        %volc_ON_no_orog = 'u-cr139'; %N.B. - not run yet so leave these
        %out for now.
        
        
        orog_str = 'orography OFF'
        emission_str = 'low altitude emissions'; 

   case '1st June runs'
 
        %um_case_PI = 'u-cw891'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF, orog ON
        um_case_PI = 'u-cx247'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano OFF, orog ON      
        um_case_PD = 'u-cx170'; run_type_DRIVER2 = 'nested ignore lat lon'; %volcano ON. orog ON    
        run_type_DRIVER = 'nested';
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
%        var_UM_DRIVER = 'Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'; 
        var_UM_DRIVER = 'LWP_sec30';
        
        %volc_OFF_no_orog = 'u-cj086';
        %volc_ON_no_orog = 'u-cr139'; %N.B. - not run yet so leave these
        %out for now.
        
        
        orog_str = 'orography OFF'
        emission_str = 'low altitude emissions'; 
        
end
        

switch SW_surf_or_TOA
    case 'TOA'
        SW_var = 'SW_up';
    otherwise
        SW_var = 'SW_down';

end

time_format_str=' UTC';

%% set the box region within which to do the cloud state analysis, timeseries, etc.
LAT_val_DRIVER2 = [30 45]; LON_val_DRIVER2 = [-60 -10]; %N Atlantic region of high negative forcing ocean only
LAT_val_DRIVER2 = [30 50]; LON_val_DRIVER2 = [-60 -10]; %N Atlantic region of high negative forcing ocean only
LAT_val_DRIVER2 = [40 50]; LON_val_DRIVER2 = [-60 -10]; %Smaller forcing region since is narrower for runs with varying SSTs
LAT_val_DRIVER2 = [33 50]; LON_val_DRIVER2 = [-60 -10]; %Extending further south for runs with varying SSTs and wind only nudging (u-ay837 and 838)

%% Set plotting region
ioverride_LAT_plots=1;
LAT_val_DRIVER_override = [8 29]; LON_val_DRIVER_override = [-180 -146]; %Whole domain Hawaii region
LAT_val_DRIVER_override = [8 25]; LON_val_DRIVER_override = [-180 -146]; %Whole width, top removed.
%LAT_val_DRIVER_override = [13 25]; LON_val_DRIVER_override = [-170 -146]; %Attempt to match MODIS image from SENSE training.
%LAT_val_DRIVER_override = [17.5 22.5]; LON_val_DRIVER_override = [-160 -152]; %Closer to Big Island
%LAT_val_DRIVER_override = [17.5 22.5]; LON_val_DRIVER_override = [-162 -152]; %Closer to Big Island - further west
%AT_val_DRIVER_override = [17.5 22.5]; LON_val_DRIVER_override = [-168 -152]; %Closer to Big Island - further west
%LAT_val_DRIVER_override = [10.5 22.5]; LON_val_DRIVER_override = [-175 -152]; %Closer to Big Island - further west


%  /home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-au652/SW_down_surf//umglaa_SW_down_surf_native res_ALL.mat    %PI
%  /home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-au536/SW_down_surf//umglaa_SW_down_surf_native res_ALL.mat     %PD

iplot_global=1; %map of SW forcing

isave_plot=0;
isave_plot_global_forcing_timser=0;
isave_plot_global_diff=0;
isave_plot_global_mean_SW=0;
isave_plot_global_lowCF=0;
isave_plot_global_pdf=0;
isave_plot_global_sw_vs_cf=0;
isave_plot_global_cont=0;
isave_plot_global_forcing=0;
isave_plot_global_midCF=0;
isave_plot_global_highCF=0;
isave_plot_global_LWP=0;
isave_plot_global_totalCF=0;

iplot_wind_arrows=0;

irestrict_domain_DRIVER=1;

clear gca

%If want to plot an outline of the nest on the global map
filename_nest = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/BC__model_level01_time08.mat';
%nest=load(filename_nest);
%eos10 may be lost... not sure what this was anyway...
irotated_pole_box=0;




scrsz=get(0,'ScreenSize');
posit=[9 60 scrsz(3) scrsz(4)];
figure('position',posit);

%% Choose times/dates here - doing all times for now
%date_str_range = 'all';

% time_round = time_ALL(it_global_diff); %this gets put on the plots
% time_round = datenum(date_str); %this gets put on the plots

% dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
% filename = [dirUM '/' run_type '_' var_UM '__ALL.mat'];
% load(filename,'time_ALL');
% load(filename,'gcm_Plon2D_UM');
%
%time_round = '01-Mar-2009';
%time_format_str = '';    

%clear time_choice; time_choice.time_range = [datenum('28-Mar-2009') datenum('29-Mar-2010')];
%clear time_choice; time_choice.time_range = [datenum('28-Mar-2009') datenum('01-Nov-2009')];
%clear time_choice; time_choice.time_range = [datenum('01-Aug-2009') datenum('31-Aug-2009')];
%clear time_choice; time_choice.time_range = [datenum('01-Dec-2009') datenum('31-Dec-2009')];
clear time_choice; time_choice.time_range = [datenum('31-Aug-2014') datenum('31-Dec-9999')];
clear time_choice; time_choice.time_range = [datenum('31-Aug-2014') datenum('04-Sep-2014 12:00')]; %currently have up to end of 0903 run for u-ba335 and 458
% Now have up to 12 UTC on 11th Sep
clear time_choice; time_choice.time_range = [datenum('31-Aug-2014') datenum('11-Sep-2014 12:00')]; %currently have up to end of 0903 run for u-ba335 and 458
clear time_choice; time_choice.time_range = [datenum('31-Aug-2014') datenum('08-Sep-2014 12:00')]; %
  %, which runs to 12UTC on 0904 (4.5 days of volcano)

  % Data starts at 12 UTC on 20th (due to spin up). Ends on 12 on 24th I
  % think.
clear time_choice; time_choice.time_range = [datenum('20-Dec-2020 12:55') datenum('24-Dec-2020 11:05')]; %reduced to not include the last timestep
%of the 4th day since for AOD, SW, etc. this is only included on the run
%starting on 24th and I did not run that far for u-ch765 (only did up to
%23rd)
%clear time_choice; time_choice.time_range = [datenum('24-Dec-2020 05:55') datenum('24-Dec-2020 11:05')]; %reduced to only not include the last timestep
%shorter period for 3d profiles to save memory

clear time_choice; time_choice.time_range = [datenum('21-Dec-2020 06:55') datenum('24-Dec-2020 11:05')]; %modified for runs where started the volcano later to reflect reality
%Currently have data for u-cj086 up to 12TUC on 25th. For u-co295 have
%until 12 UTC on 26th Dec. So could extend out to 12UTC on 25th if needed.

%clear time_choice; time_choice.time_range = [datenum('20-Dec-2020') datenum('23-Dec-2020 13:00')]; %temp shorter timerseries

%Specific days :-
%clear time_choice; time_choice.time_specific = datenum('06-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('07-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('10-Aug-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Dec-2009'); time_choice.find_nearest=1;


clear time_choice; time_choice.time_range = [datenum('01-Jun-2018 00:00') datenum('15-Jun-2018 00:05')]; %modified for runs where started the volcano later to reflect reality


%% Get time data from the specified file
    var_UM = var_UM_DRIVER;
    
    switch run_type_DRIVER
        case 'nested'
            um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = 'nested'; icoarse=0; ivar_dir=1; %nested UKCA
            dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
            dirUM_PD = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case_PD];
            
        case 'nested ignore lat lon'
            um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = 'nested ignore lat lon'; icoarse=0; ivar_dir=1; %nested UKCA
            dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];            

        case 'global'
            % PI run
            %    um_case='u-au652'; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
            um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
            dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];

    end

    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
    dat_global_PD = UM_load_merged_netCDF(dirUM_PD,var_UM,run_type,load_type,pole_lat,pole_lon);

    

array_in=[]; %just test to get the indices for now.
dim=NaN; %don't need dim if just getting the indices
[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);
%Also do for the PD run in case the times are mismatched (requiring
%different indices).
[out, time_out_PD, time_inds_PD, dtime_match_PD] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);
if length(time_inds)==1
    time_round = datestr(dat_global.time_ALL(time_inds(1)));
else
    time_round = [datestr(dat_global.time_ALL(time_inds(1))) ' to ' datestr(dat_global.time_ALL(time_inds(end)))];
end

[Y_UM,M_UM,D_UM,HH_UM,MM_UM,SS_UM] = datevec(time_out);


    
    gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
    gcm_Plat2D_edges_UM = dat_global.gcm_Plat2D_edges_UM;

    gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
    %    if maxALL(gcm_Plon2D_UM) < 181
    if irestrict_domain_DRIVER==0
        i180 = find(gcm_Plon2D_UM<0);
        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180)+360;
    end
    
    switch run_type
        case 'nested ignore lat lon'
        otherwise
            %[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
            [gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
            %    gcm_Plon2D_edges_UM = dat_global.gcm_Plon2D_edges_UM;
    end

% %    nT = length(eval(['dat_global.' var_UM '_ALL;']));
%     nT = length(time_inds);
%     dat_modis_PI_ALL = NaN * ones([size(gcm_Plat2D_UM,1) size(gcm_Plat2D_UM,2) nT]);    
%     it=0;
%     for it_global_diff=time_inds
%         it=it+1;
%         dat_modis_PI_ALL(:,:,it) = eval(['dat_global.' var_UM '_ALL{it_global_diff};']);        
%     end

save_time_lat_file = [dirUM '/save_time_lat.mat'];
save(save_time_lat_file,'-V7.3');
    
 
%% load SW_down_TOA
iload_SW_down_TOA=0;
if iload_SW_down_TOA==1
%TOA SW downwards (incoming) - N.B. - this is offset by 20 mins from the
%other SW diags (e.g. SW_down_surf) and so needs to have its own time
%indicies caculated for it. This gives the same number of time datapoints as before. Presumably this is because it is output on
%radiation timesteps.
var_UM = 'SW_down_TOA';
 
um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);

array_in=[]; %just test to get the indices for now.
dim=NaN; %don't need dim if just getting the indices
[out, time_out_SW_in, time_inds_SW_in, dtime_match_SW_in] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

[SW_in_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_SW_in,load_type,gcm_Plat2D_UM);

um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
[SW_in_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_SW_in,load_type,gcm_Plat2D_UM);

    
end   



%% Load satellite data, etc.

%CALIOSP CF data - load using script :- read_calipso_monthly_night_IPSL.m for average values
    %first. This puts the data into these fields :- 
    % cllcalipso_monthly_AVERAGE, clmcalipso_monthly_AVERAGE, clhcalipso_monthly_AVERAGE
    % Is average of day and night values for now (have separate ones too).
if ical_data==1
    read_calipso_monthly_night_IPSL
    [out, time_out_CAL, time_inds_CAL, dtime_match] = get_time_range_of_array(array_in,gcm_time_matlab_CALIPSO_monthly,time_choice,dim);
end

if iceres_data==1
    ceres_file = '/home/disk/eos15/d.grosvenor/eos8/CERES/ACSIS/CERES_EBAF-TOA_Ed4.0_Subset_200903-201004.nc';
        %Monthly averages from March 2009 to April 2010.
    nc_ceres = netcdf(ceres_file);
    lon_ceres = nc_ceres{'lon'}(:);
    lat_ceres = nc_ceres{'lat'}(:);    
    time_ceres = nc_ceres{'time'}(:); %days since 2000-03-01
    SW_TOA_ceres = nc_ceres{'toa_sw_all_mon'}(:);
    
    time_matlab_CERES = datenum('01-Mar-2000') + time_ceres;
    [out, time_out_CERES, time_inds_CERES, dtime_match] = get_time_range_of_array(array_in,time_matlab_CERES,time_choice,dim);
    
    i180=find(lon_ceres>180);
    lon_ceres(i180) = lon_ceres(i180)-360;
    [gcm_Plon2D_CERES, gcm_Plat2D_CERES] = meshgrid(lon_ceres,lat_ceres);
    
end


if amsre_data==1
    %Global AMSRE data saved for 2005 to 2011 (1st Dec 2005 to 31st Jan
    %2011)
%    amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_2005_to_2011_20130812T220913.mat';
    amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_Aqua_2005_to_2011_20180814T134447.mat' %inc. rain rates, wind speeds, etc.
    if ~exist('amsre_dat')
        amsre_dat = load(amsre_matfile,'lwp_amsre','day_amsre','month_amsre','year_amsre','rain_amsre','sst_amsre');
    end
    amsre_matlab_time = datenum(amsre_dat.year_amsre,amsre_dat.month_amsre,amsre_dat.day_amsre);
    [out, time_out_AMSRE, time_inds_AMSRE, dtime_match] = get_time_range_of_array(array_in,amsre_matlab_time,time_choice,dim);
    
    
    %convert rain rate and LWP into LWP + RWP
end


%%
error('*** ENDING here - No errors, just the end of the loading stuff! ***');

%% SW direct and indirect forcing 
isave_plot_global_diff=0;
isave_plot_global_forcing=0;



% Will just consider clouds in the PI run for now
% consider coarse graining to help make sure met is the same?

%Clean means no direct effects (background clean aerosol when out of
%clouds)
% Direct effect is the total (indirect + direct) minus SW_clean (just
% indirect effect, but no direct), i.e. direct = indirect+direct - indirect

 var_UM = [SW_var '_' SW_surf_or_TOA];
 
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [SW_down_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [SW_down_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    
% Do some plots of this       
    dat_modis_PI = meanNoNan(SW_down_PI_ALL,3);
    
    %restrict to the region of interest
    if i_select_region==1
        LAT_val = LAT_val_DRIVER2;
        LON_val = LON_val_DRIVER2;

        [iregion_lin,iregion_lin_edges,SW_down_region_PI]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,SW_down_PI_ALL(:,:,1:end));
        SW_down_region_PI_timser = meanNoNan(meanNoNan(SW_down_region_PI,1),1);
    end
            
    
    dat_modis_PD = meanNoNan(SW_down_PD_ALL,3);
    
    %restrict to the region of interest
    if i_select_region==1
        LAT_val = LAT_val_DRIVER2;
        LON_val = LON_val_DRIVER2;

        [iregion_lin,iregion_lin_edges,SW_down_region_PD]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,SW_down_PD_ALL(:,:,1:end));
        SW_down_region_PD_timser = meanNoNan(meanNoNan(SW_down_region_PD,1),1);
    end


    dat_modis = dat_modis_PD - dat_modis_PI;
    
    %run plotting script
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
    caxis([-25 25]);
    titlenam_driver='\DeltaSW_{TOA}';
    title(titlenam_driver);
    increase_font_size_map_figures

% Save
%savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/model_vs_MODIS_Nd_plots/' modis_day_of_year_str '_global_' titlenam_driver];
%savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
savename=['/home/disk/eos1/d.grosvenor/modis_work/Iceland/UKCA/difference_plots/' run_set '_' titlenam_driver];

if isave_plot_global_diff==1
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end
    

 var_UM = [SW_var '_clean_' SW_surf_or_TOA]; 
 
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [SW_nodirect_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [SW_nodirect_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);    
    
    aerosol_effect_PI_ALL = SW_down_PI_ALL - SW_nodirect_PI_ALL;    
    aerosol_effect_PD_ALL = SW_down_PD_ALL - SW_nodirect_PD_ALL;
    direct_ALL = aerosol_effect_PD_ALL - aerosol_effect_PI_ALL;
    
    
% For indirect effect we do clean minus clear_clean.
% I.e., no aerosol effect with cloud effect - no aerosol effect without cloud effect
    %no direct or indirect effects
    
 var_UM = [SW_var '_clean_clear_'  SW_surf_or_TOA]; 
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [SW_clean_clear_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [SW_clean_clear_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    
    cloud_effect_PI_ALL = SW_nodirect_PI_ALL - SW_clean_clear_PI_ALL;    
    cloud_effect_PD_ALL = SW_nodirect_PD_ALL - SW_clean_clear_PD_ALL;
    indirect_ALL = cloud_effect_PD_ALL - cloud_effect_PI_ALL;
    
    switch SW_surf_or_TOA
        case 'TOA'
            indirect_ALL = - indirect_ALL;
    end
     
        
%    dat_modis = meanNoNan(low_CF_PD_ALL,3) - meanNoNan(low_CF_PI_ALL,3);
    
    var_UM = 'Direct forcing (W m^{-2})';
    dat_modis = meanNoNan(direct_ALL,3);    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-10 10]);
    
    % Save
    %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/model_vs_MODIS_Nd_plots/' modis_day_of_year_str '_global_' titlenam_driver];
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];

    if isave_plot_global_forcing==1
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        close(gcf);
    end
    
%% Plot indirect
isave_plot_global_diff=0;

    var_UM = 'Indirect forcing (W m^{-2})';
    dat_modis = meanNoNan(indirect_ALL,3);  
    
    %Save the forcing for Dan Jones
    time_range = 24*(datenum(time_choice.time_range - datenum('01-Jan-1970'))); %Hrs since 1st Jan 1970
    indirect_forcing_time_mean = dat_modis;
    lat2D = gcm_Plat2D_UM;
    lon2D = gcm_Plon2D_UM;
    
    if isave_Dan_Jones==1
        save_file_Dan_Jones = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/indirect_forcing_timemean.mat'];
        save(save_file_Dan_Jones,'indirect_forcing_time_mean','lat2D','lon2D','time_range','-V7.3');
        mat2nc_Dan(save_file_Dan_Jones,[save_file_Dan_Jones '.nc']);
    end
    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    increase_font_size_map_figures
     caxis([-15 15]);
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    
    % Save  
    if isave_plot_global_forcing==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/Iceland/' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end
    
%Also the CRE (SW_down due to clouds at surface)    

    var_UM = 'Indirect cloud effect PI (W m^{-2})';
    dat_modis = meanNoNan(cloud_effect_PI_ALL,3);    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-120 0]);
    
    % Save  
    if isave_plot_global_forcing==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        close(gcf);
    end
    
    %PD    
    
    var_UM = 'Indirect cloud effect PD (W m^{-2})';
    dat_modis = meanNoNan(cloud_effect_PD_ALL,3);    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-120 0]);
    
    % Save  
    if isave_plot_global_forcing==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        close(gcf);
    end


%% Make a timeseries
isave_plot_global_forcing_timser=0;

    %(restricting to the region of interest)
    [indirect_ALL_timser] = UM_make_regional_timeseries(indirect_ALL,nT,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM);
    
    figure
    time = dat_global.time_ALL(time_inds);
    plot(time,indirect_ALL_timser,'linewidth',3);
    datetick('x','mmm');
    
    %legend(leg_strs,'location','northwest');
    xlabel('Time');
    %ylabel('SW surface forcing (W m^{-2})');
    ylabel('SW surface indirect forcing (W m^{-2})');
    fontsize_figure(gcf,gca,18);
%    set(gca,'xlim',[-0.05 1.05]);

% Save
    if isave_plot_global_forcing_timser==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_indirect_SW_timser'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        close(gcf);
    end

%% Load SW_up_TOA and 
%% calculate change in SW out TOA and SW out TOA vs CERES
isave_plot_CERES=1;

    clear gca
    var_UM = 'SW_up_TOA';
 
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [SW_TOA_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = 'nested ignore lat lon';
    %run_type_DRIVER; 
    icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [SW_TOA_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    dat_modis = meanNoNan(SW_TOA_PD_ALL,3) - meanNoNan(SW_TOA_PI_ALL,3); var_UM = 'Change in SW out TOA (PD minus PI)';
    subtitle_str = var_UM;
%    dat_modis = meanNoNan(SW_TOA_PD_ALL,3);
%    dat_modis = meanNoNan(SW_TOA_PI_ALL,3);
   

  isave=1;
    if isave==1
        CF_str_var = 'SW_TOA';
        
        save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PI_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PD_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    end
    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%    caxis([-0.1 0.1]);   


if iceres_data==1   
    ceres_SW_TOA_time_period_mean = meanNoNan(SW_TOA_ceres(time_inds_CERES,:,:),1);
    UM_ACSIS_global_SUBPLOT_commands_CERES
end


%% Plot cloud fractions - low cloud
isave_plot_global_lowCF=0;
isave_plot_CF = 0;

clear gca

switch cloud_input
    case 'UM'
        var_UM = 'low_cloud_amount';
        icosp_mask=0;
    case 'CALIPSO'
        var_UM = 'calipso_low_cloud_amount';
        var_UM_mask = [var_UM '_mask'];
        icosp_mask=1;
    case 'MODIS'
        var_UM = 'modis_low_cloud_amount';
        var_UM_mask = ['modis_misr_issp_cloud_amount_mask'];
        icosp_mask=1;
end
 
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [low_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    low_CF_PD_ALL(low_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??  
    
    if icosp_mask==1        
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [low_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        low_CF_PD_ALL(low_CF_mask==0) = NaN;
    end
    
    %var_UM = 'low_cloud_amount';
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [low_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    low_CF_PI_ALL(low_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??
    
    low_CF_av_ALL = (low_CF_PI_ALL + low_CF_PD_ALL) /2;
    
    
    isave=0;
    if isave==1
        CF_str_var = 'low_CF';
        
        save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PI_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PD_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    end
        
    
    %% dlowCF - Restrict region to threshold DeltaNd - filtered for each time separately
    thresh_Nd = 600; %per cc
    dNd = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
    inan = find(dNd < thresh_Nd); 
    
    dCF = low_CF_PD_ALL - low_CF_PI_ALL;    
    dCF_plume = dCF;
    dCF_plume(inan)=NaN;
    dom_mean_dCF_plume = meanNoNan(meanNoNan(dCF_plume,1),1);           
    
    figure
    plot(time,dom_mean_dCF_plume,'linewidth',3);
    datetick('x','dd');
    
    %legend(leg_strs,'location','northwest');
    xlabel('Time');
    %ylabel('SW surface forcing (W m^{-2})');
    ylabel('\DeltalowCF (g m^{-2})');
    set(gca,'ylim',[-0.15 0.15]);
    fontsize_figure(gcf,gca,18);
    grid on
    title(['\DeltaN_d thresh=' num2str(thresh_Nd) ' cm^{-3}']);
    
 %% dlowCF - Inside vs outside plume for just volcano ON
    %Restrict region to threshold DeltaNd - filtered for each time separately
    
    thresh_Nd = 200; %per cc
    dNd = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
    inan = find(dNd < thresh_Nd); 
        
    CF_plume = low_CF_PD_ALL;
    CF_plume(inan)=NaN;
    dom_mean_CF_plume = meanNoNan(meanNoNan(CF_plume,1),1);           
    
    CF_non_plume = NaN*ones(size(low_CF_PD_ALL));
    CF_non_plume(inan)=low_CF_PD_ALL(inan);
    dom_mean_CF_non_plume = meanNoNan(meanNoNan(CF_non_plume,1),1);    
    
    figure
    plot(time,dom_mean_CF_plume-dom_mean_CF_non_plume,'linewidth',3);
    datetick('x','dd');
    
    %legend(leg_strs,'location','northwest');
    xlabel('Time');
    %ylabel('SW surface forcing (W m^{-2})');
    ylabel('\DeltalowCF, plume minus non-plume');
    set(gca,'ylim',[-0.1 0.15]);
    fontsize_figure(gcf,gca,18);
    grid on
    title(['\DeltaN_d thresh=' num2str(thresh_Nd) ' cm^{-3}']);    
    
%% 
    dat_modis = meanNoNan(low_CF_PD_ALL,3) - meanNoNan(low_CF_PI_ALL,3); var_UM = 'Change in low CF (PD minus PI)';
%    dat_modis = meanNoNan(low_CF_PD_ALL,3);
%    dat_modis = meanNoNan(low_CF_PI_ALL,3);
   
    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-0.1 0.1]);
    
    
    % Save
    if isave_plot_global_lowCF==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        close(gcf);
    end
    
    %Plot CALIPSO data for same period
    %Run script :- read_calipso_monthly_night_IPSL.m for average values
    %first. This puts the data into these fields :- 
    % cllcalipso_monthly_AVERAGE, clmcalipso_monthly_AVERAGE, clhcalipso_monthly_AVERAGE
    % Is average of day and night values for now (have separate ones too).

if ical_data==1    
    var_UM = 'CALIPSO low cloud fraction';
    dat_modis = meanNoNan(cllcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
    figure
    %run plotting script
    cloud_alt = 'Low altitude cloud';
    CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 1]);  
    
        
    %Run subplotting script for UM, CALIPSO and bias
    UM_ACSIS_global_SUBPLOT_commands_cloud_eval
    
end
    
    
%% mid cloud
isave_plot_global_midCF=0;
isave_plot_CF = 0;

clear gca    

switch cloud_input
    case 'UM'
        var_UM = 'mid_cloud_amount';
        icosp_mask=0;
    case 'CALIPSO'
        var_UM = 'calipso_mid_cloud_amount';
        var_UM_mask = [var_UM '_mask'];
        icosp_mask=1;
    case 'MODIS'
        var_UM = 'modis_mid_cloud_amount';
        var_UM_mask = ['modis_misr_issp_cloud_amount_mask'];
        icosp_mask=1;
end

    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [mid_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    mid_CF_PD_ALL(mid_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [mid_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        mid_CF_PD_ALL(mid_CF_mask==0) = NaN;
    end
    
%Don't have COSP output for PI yet
    var_UM = 'mid_cloud_amount';
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);     
    [mid_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    mid_CF_PI_ALL(mid_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??   
    
    mid_CF_av_ALL = (mid_CF_PI_ALL + mid_CF_PD_ALL) /2;
    
    
    isave=1;
    if isave==1
        CF_str_var = 'mid_CF';
        
        save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PI_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PD_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    end
    

    %% Plot cloud fractions - mid cloud
    %dat_modis = meanNoNan(mid_CF_PI_ALL,3);
    %dat_modis = meanNoNan(mid_CF_PD_ALL,3); 
    dat_modis = meanNoNan(mid_CF_PD_ALL,3) - meanNoNan(mid_CF_PI_ALL,3);
    
    var_UM = 'Change in mid CF (PD minus PI)';

    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-0.1 0.1]);

    
    % Save
    if isave_plot_global_midCF==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        close(gcf);
    end
    
    
    if ical_data==1

        dat_modis = meanNoNan(clmcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
        figure
        %run plotting script
        cloud_alt = 'Mid altitude cloud';
        CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 1]);

        %Run subplotting script for UM, CALIPSO and bias
        UM_ACSIS_global_SUBPLOT_commands_cloud_eval

    end
    
    
%% High cloud        
isave_plot_global_highCF=0;
isave_plot_CF = 0;

clear gca

    switch cloud_input
        case 'UM'
            var_UM = 'high_cloud_amount';
            icosp_mask=0;
        case 'CALIPSO'
            var_UM = 'calipso_high_cloud_amount';
            var_UM_mask = [var_UM '_mask'];
            icosp_mask=1;
        case 'MODIS'
            var_UM = 'modis_high_cloud_amount';
            var_UM_mask = ['modis_misr_issp_cloud_amount_mask'];
            icosp_mask=1;
    end
    
 
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [high_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    high_CF_PD_ALL(high_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [high_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        high_CF_PD_ALL(high_CF_mask==0) = NaN;
    end
    
    var_UM = 'high_cloud_amount';
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon); 
    [high_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    high_CF_PI_ALL(high_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??   
    
    high_CF_av_ALL = (high_CF_PI_ALL + high_CF_PD_ALL) /2;
    
    isave=1;
    if isave==1
        CF_str_var = 'high_CF';
        
        save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PI_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PD_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    end

       
    
    %dat_modis = meanNoNan(high_CF_PI_ALL,3);
    %dat_modis = meanNoNan(high_CF_PD_ALL,3);    
    dat_modis = meanNoNan(high_CF_PD_ALL,3) - meanNoNan(high_CF_PI_ALL,3); var_UM = 'Change in high CF (PD minus PI)';
    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-0.1 0.1]);
        
    
    % Save


    if isave_plot_global_highCF==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        close(gcf);
    end
    
    if ical_data==1
                
        dat_modis = meanNoNan(clhcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
        figure
        %run plotting script
        cloud_alt = 'High altitude cloud';
        CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 1]);


        %Run subplotting script for UM, CALIPSO and bias
        UM_ACSIS_global_SUBPLOT_commands_cloud_eval

    end

    
%% Total cloud   
isave_plot_global_totalCF=0;

icalc_total=1;

    switch cloud_input
        case 'UM'
            var_UM = 'NOT DEFINED?'; %Can get this from stash 9-217 (is added in Python script now).
            icosp_mask=0;                                    
            icalc_total=0;
        case 'CALIPSO'
            var_UM = 'calipso_total_cloud_amount';
            var_UM_mask = [var_UM '_mask'];
            icosp_mask=1;
        case 'MODIS'
            var_UM = 'modis_total_cloud_amount';
            var_UM_mask = ['modis_misr_issp_cloud_amount_mask'];
            icosp_mask=1;
    end
    
    if icalc_total==1

        um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
        [total_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??

        if icosp_mask==1
            dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
            [total_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
            total_CF_PD_ALL(total_CF_mask==0) = NaN;
        end

        %     var_UM = 'total_cloud_amount';
        %     um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        %     dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        %     dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
        %     [total_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        %     total_CF_PI_ALL(total_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??
        %
        %
        
        %quick plot of low+mid+high for comparing cyclone positions - PD
        var_UM = 'PD low + mid + high cloud fraction'
        dat_modis = meanNoNan(low_CF_PD_ALL,3) + meanNoNan(mid_CF_PD_ALL,3) + meanNoNan(high_CF_PD_ALL,3);

        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 2]);
                        
        
        %quick plot of low+mid+high for comparing cyclone positions - PI
        var_UM = 'PI low + mid + high cloud fraction'
        dat_modis = meanNoNan(low_CF_PI_ALL,3) + meanNoNan(mid_CF_PI_ALL,3) + meanNoNan(high_CF_PI_ALL,3);

        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 2]);

        % Save
        if isave_plot_global_totalCF==1
            savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
            clear opts
            %        opts.iplot_png=1;
            opts.iplot_eps=1;
            saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
            close(gcf);
        end

        if ical_data==1

            %     dat_modis = meanNoNan(clhcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
            %     figure
            %     %run plotting script
            %     CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
            %     lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
            %     caxis([0 1]);


            %Run subplotting script for UM, CALIPSO and bias
            cloud_alt = 'Total cloud';
            UM_ACSIS_global_SUBPLOT_commands_cloud_eval

        end

    end
    
%% Droplet concs    
isave_plot_global_Nd = 0;
isave_plot_Nd_bias=1;
iload_modis_dat=0;    

if iload_modis_dat==1
    % saved some global MODIS (mockL3) data for the period [datenum('28-Mar-2009') datenum('29-Mar-2010')]; here :-
    %CF>80
    %modis_Nd_file = '/home/disk/eos10/d.grosvenor/saved_MODIS_data/NAtlantic/global_mockL3_CF80_SZA65_2009-2010_ACSIS.mat';
    modis_Nd_file = '/home/disk/eos10/d.grosvenor/saved_MODIS_data/Iceland/global_mockL3_CF80_SZA65_30-Aug-10-Sep-2014_ACSIS.mat';
    %CF>0
    % modis_Nd_file = '/home/disk/eos10/d.grosvenor/saved_MODIS_data/NAtlantic/global_mockL3_CF0_SZA65_2009-2010_ACSIS.mat';
    % Load on first use :-
    if ~exist('Cloud_Fraction_Liquid')        
        modis_loaded = load(modis_Nd_file,'CTH','W_time3','Cloud_Fraction_Liquid','Droplet_Number_Concentration_37',...
         'gcm_Plat2D_AMSRE','gcm_Plon2D_AMSRE','modisyear_timeseries3','daynum_timeseries3');
    end
    
end

 Nd_type = 'lwc weighted';
 Nd_type = 'lwc in-cloud weighted'; 
 Nd_type = 'CASIM div CF lwc in-cloud weighted z3000'; 
 Nd_type = 'CASIM div CF lwc in-cloud weighted domain top'; 
 
    switch Nd_type
        case 'lwc weighted'
            var_UM2 = ['Nd_lwc_weighted_' Nd_var_str];
            var_UM = [var_UM2 Nd_var_str2];
        case 'lwc in-cloud weighted'
            var_UM2 = ['Nd_lwc_in_cloud_weighted_' Nd_var_str];     
            var_UM = [var_UM2 Nd_var_str2];
        case 'CASIM div CF lwc in-cloud weighted z3000'
            var_UM = ['Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_z3000'];
            var_UM2 = var_UM;
        case 'CASIM div CF lwc in-cloud weighted domain top'
            var_UM = ['Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'];
            var_UM2 = var_UM;
            
    end
    
    icosp_mask=0;
    
        
 
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon,var_UM2); 
    clear opts; opts.time_out_check = time_out;
    [time_inds_var,time_out_var] = UM_get_time_inds_from_dat_global(dat_global,time_choice,opts);   
    [Nd_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_var,load_type,gcm_Plat2D_UM);  %per m3
%    total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [total_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        total_CF_PD_ALL(total_CF_mask==0) = NaN;
    end
    

    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon,var_UM2);  
    clear opts; opts.time_out_check = time_out;
    [time_inds_var,time_out_var] = UM_get_time_inds_from_dat_global(dat_global,time_choice,opts);  
    [Nd_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_var,load_type,gcm_Plat2D_UM);  %per m3
%    total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
 
    
    

    
%     dat_modis = meanNoNan(clhcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
%     figure
%     %run plotting script
%     CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%     lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%     caxis([0 1]);
    
if iload_modis_dat==1    
    %Run subplotting script for UM, MODIS and bias
    UM_ACSIS_global_SUBPLOT_commands_Nd
end

    save_file = [UM_base_dir um_case_PI '/Nd.mat'];
    Hawaii_aod_coarse_grain_save(save_file,Nd_PI_ALL,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    
    save_file = [UM_base_dir um_case_PD '/Nd.mat'];
    Hawaii_aod_coarse_grain_save(save_file,Nd_PD_ALL,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
%Change in Nd from PI to PD

    Nd_PI_map = meanNoNan(Nd_PI_ALL,3);
    Nd_PD_map = meanNoNan(Nd_PD_ALL,3);
    
   
    
    
    dat_modis = (Nd_PD_map - Nd_PI_map)/1e6; var_UM = 'Change in Nd (PD minus PI; cm^{-3})';  
    
    iplot_Nd=1;
    if iplot_Nd==1
        %run plotting script
        figure
        subtitle_str = var_UM;
        ilabel_colorbar = 1; col_bar_lab_str = '(cm ^{-3})';        
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));        
        caxis([0 300]);
        clear ilabel_colorbar        
    end
    
    it=96; dat_modis = (Nd_PD_ALL(:,:,it) - Nd_PI_ALL(:,:,it))/1e6; subtitle_str = ['Change in Nd (PD minus PI; cm^{-3}) it=' num2str(it)];  
    var_UM='';
    iplot_Nd=1;
    if iplot_Nd==1
        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-100 100]);
        
    end
    
    
    delta_Nd_map = dat_modis;
    
    % Save

    if isave_plot_global_Nd==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end
    
    dat_modis = 100 * (Nd_PD_map - Nd_PI_map) ./ Nd_PI_map ; var_UM = 'Percentage change in Nd (PD minus PI)';    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-150 150]);  
    
    
    if isave_plot_global_Nd==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end
    
    it=96; dat_modis = Nd_PI_ALL(:,:,it)/1e6; subtitle_str = ['Volcano OFF Nd (cm^{-3}) it=' num2str(it)];  
    var_UM='';
    iplot_Nd=1;
    if iplot_Nd==1
        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 100]);
        
    end

    
%% UKCA Nd

 Nd_type = 'lwc weighted'; %UKCA Nd
 Nd_type = 'lwc in-cloud weighted'; %UKCA Nd
 %Nd_type = 'CASIM div CF lwc in-cloud weighted z3000'; 
 %Nd_type = 'CASIM div CF lwc in-cloud weighted domain top'; 
 
    switch Nd_type
        case 'lwc weighted' %UKCA Nd - calculated diagnostically and no need to div by CF.
            var_UM2 = ['Nd_lwc_weighted_' Nd_var_str];
            var_UM = [var_UM2 Nd_var_str2];
        case 'lwc in-cloud weighted' %UKCA Nd
            var_UM2 = ['Nd_lwc_in_cloud_weighted_' Nd_var_str];     
            var_UM = [var_UM2 Nd_var_str2];
        case 'CASIM div CF lwc in-cloud weighted z3000'
            var_UM = ['Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_z3000'];
            var_UM2 = var_UM;
        case 'CASIM div CF lwc in-cloud weighted domain top'
            var_UM = ['Nd_div_CF_lwc_in_cloud_weighted_CASIM_total_column_to_zdomain_top'];
            var_UM2 = var_UM;
            
    end
    
    icosp_mask=0;
             
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon,var_UM2);   
    [Nd_UKCA_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
%    total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [total_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        total_CF_PD_ALL(total_CF_mask==0) = NaN;
    end
    

    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon,var_UM2);   
    [Nd_UKCA_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
%    total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    
    
%% Load and plot LWP
isave_plot_global_LWP=0;
isave_plot_LWP_bias=0;
iplot_amsre=0;

clear gca

 var_UM_LWP = 'LWP'; %calculated from Python script using QCL after timestep(0-254)
 var_UM_LWP = 'LWP_sec30'; %from section 30-405

    var_UM = var_UM_LWP;
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [LWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    var_UM = var_UM_LWP; 
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [LWP_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);   

    if exist('volc_OFF_no_orog')
	    %No orog runs
	    var_UM = var_UM_LWP;
	    um_case=volc_OFF_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
	    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
	    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
	    [LWP_PI_ALL2,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
	    
	    var_UM = var_UM_LWP; 
	    um_case=volc_ON_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER2; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
	    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
	    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
	    [LWP_PD_ALL2,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  
	    
    else
        LWP_PI_ALL2 = NaN;
        LWP_PD_ALL2 = NaN;        
    end
    
    
    switch var_UM
        case {'LWP_sec30','LWP_correct_dz'}
           LWP_PD_ALL =  LWP_PD_ALL*1e3; %convert to g/m2
           LWP_PI_ALL =  LWP_PI_ALL*1e3; %convert to g/m2  
           LWP_PD_ALL2 =  LWP_PD_ALL2*1e3; %convert to g/m2
           LWP_PI_ALL2 =  LWP_PI_ALL2*1e3; %convert to g/m2  
    end
    
    mean_LWP_PI = meanNoNan(LWP_PI_ALL,3);
    mean_LWP_PD = meanNoNan(LWP_PD_ALL,3); 
    mean_LWP_PI2 = meanNoNan(LWP_PI_ALL2,3);
    mean_LWP_PD2 = meanNoNan(LWP_PD_ALL2,3); 
    
    isave=0;
    if isave==1
        CF_str_var = 'LWP';
        
        save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval(['LWP_PI_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval(['LWP_PD_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    end

    
    icoarse_grain=1;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = 'Time mean \DeltaLWP (g m^{-2})';
    var_UM = '';
    dat_modis = mean_LWP_PD - mean_LWP_PI;
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]);
    icoarse_grain=0;
   
    icoarse_grain=1;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = 'Time mean % change in all-sky LWP (g m^{-2})';
    var_UM = '';
    dat_modis = 100*(mean_LWP_PD - mean_LWP_PI)./mean_LWP_PI;
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-30 30]);
    icoarse_grain=0;

        isave_plot=1;
        if isave_plot==1
                savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_LWP_' titlenam_driver '_vs_satellite'];
                clear opts
                        opts.iplot_png=1;
                opts.iplot_eps=1;
                saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
                %    close(gcf);
        end


 
    icoarse_grain=1;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = 'Time mean \DeltaLWP (g m^{-2}) \tau>3, LWP>3';
    var_UM = '';
    LWP_PI = LWP_PI_ALL;
    %inan = find(tau_PI_ALL<3);
    inan = find(tau_PI_ALL<3 | LWP_PI_ALL<3);
    LWP_PI(inan) = NaN;    
    LWP_PD = LWP_PD_ALL;
    %inan = find(tau_PD_ALL<3);
    inan = find(tau_PD_ALL<3 | LWP_PD_ALL<3);
    LWP_PD(inan) = NaN;    
    dat_modis = meanNoNan(LWP_PD,3) - meanNoNan(LWP_PI,3);
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]);
    icoarse_grain=0;
    
    it=96;
    dat_modis = LWP_PD_ALL(:,:,it) - LWP_PI_ALL(:,:,it);
    icoarse_grain=1;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['Time mean \DeltaLWP (g m^{-2}), it=' num2str(it)];
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-10 10]);
    icoarse_grain=0;
    
    it=96;
    dat_modis = LWP_PI_ALL(:,:,it);
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['LWP Volcano OFF (g m^{-2}), it=' num2str(it)];
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 400]);
    icoarse_grain=0;
    
    it=23;
    it=22;
    it=33;
    dat_modis = LWP_PD_ALL(:,:,it);
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['LWP Volcano ON (g m^{-2}), ' datestr(time_out(it))];
    ilabel_colorbar = 1; col_bar_lab_str = '(g m^{-2})';
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 200]);
    icoarse_grain=0;
    
    
    if exist('W0_orig')
        it=33;
        dat_modis = f1_orig(:,:,it);
        icoarse_grain=0;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = ['CF_{subgrid} Volcano ON (g m^{-2}), ' datestr(time_out(it))];
        ilabel_colorbar = 1; col_bar_lab_str = '';
        var_UM = '';
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 200]);
        icoarse_grain=0;
        
        icoarse_grain=0;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = 'Time mean \DeltaLWP_{ic subgrid} (g m^{-2})';
        var_UM = '';
        dat_modis = meanNoNan(W1_orig,3) - meanNoNan(W0_orig,3);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '(g m^{-2})';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global        
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-50 50]);
        icoarse_grain=0;
        
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = 'Time mean % \DeltaLWP_{ic subgrid} (g m^{-2})';
        var_UM = '';
        dat_modis = 100*(meanNoNan(W1_orig,3) - meanNoNan(W0_orig,3)) ./ meanNoNan(W0_orig,3);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '(%)';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-50 50]);
        icoarse_grain=0;
        
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = 'Time mean \DeltaCF_{subgrid} (g m^{-2})';
        var_UM = '';
        dat_modis = meanNoNan(f1_orig,3) - meanNoNan(f0_orig,3);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-0.5 0.5]);
        icoarse_grain=0;
        
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = 'Time mean % \DeltaCF_{subgrid} (%)';
        var_UM = '';
        dat_modis = 100*(meanNoNan(f1_orig,3) - meanNoNan(f0_orig,3)) ./ meanNoNan(f0_orig,3);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '(%)';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-20 20]);
        icoarse_grain=0;
        
                icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = 'Time mean % \DeltaCF_{subgrid} (%)';
        var_UM = '';
        dat_modis = 100*(meanNoNan(f1_orig,3) - meanNoNan(f0_orig,3)) ./ meanNoNan(f0_orig,3);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '(%)';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-20 20]);
        icoarse_grain=0;
        
%        total_max_random_cloud_amount
        
        icoarse_grain=0;
        M_coarse_grain=10; N_coarse_grain=10;
        subtitle_str = 'Time mean CF_{subgrid}';
        var_UM = '';
        dat_modis = meanNoNan(f1_orig,3);
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 1]);
        icoarse_grain=0;
        
        
    end
        


    
%% time mean LWP plots (not deltas) - orograpy and no orog

%Volcano ON.
    dat_modis = mean_LWP_PD;
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['LWP Volcano ON (g m^{-2}), time mean'];
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 200]);
    icoarse_grain=0;
    
    dat_modis = mean_LWP_PD2;
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['LWP Volcano ON, orog OFF (g m^{-2}), time mean'];
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 200]);
    icoarse_grain=0;
    
    
%Volcano OFF.
    dat_modis = mean_LWP_PI;
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['LWP Volcano OFF (g m^{-2}), time mean'];
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 200]);
    icoarse_grain=0;
    
    dat_modis = mean_LWP_PI2;
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['LWP Volcano OFF, orog OFF (g m^{-2}), time mean'];
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 200]);
    icoarse_grain=0;    
    
    
%% more LWP plots
    
    %save the timemean values to disk
    %save_time_lat_file = [dirUM '/save_time_lat.mat'];
    save(save_time_lat_file,'mean_LWP_PI','mean_LWP_PD','-APPEND','-V7.3');
    
    
%Plot timeseries of the mean difference
    dom_mean_LWP_PI = meanNoNan(meanNoNan(LWP_PI_ALL,1),1);
    dom_mean_LWP_PD = meanNoNan(meanNoNan(LWP_PD_ALL,1),1);
    
    save(save_time_lat_file,'dom_mean_LWP_PI','dom_mean_LWP_PD','-APPEND','-V7.3');
    
    time_mean_dLWP = meanNoNan(dom_mean_LWP_PD-dom_mean_LWP_PI,1);
    
%    [indirect_ALL_timser] = UM_make_regional_timeseries(indirect_ALL,nT,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM);
        
    time = dat_global.time_ALL(time_inds);
    
    figure
    plot(time,dom_mean_LWP_PD-dom_mean_LWP_PI,'linewidth',3);
    datetick('x','dd');
    
    %legend(leg_strs,'location','northwest');
    xlabel('Time');
    %ylabel('SW surface forcing (W m^{-2})');
    ylabel('\DeltaLWP (g m^{-2})');
    set(gca,'ylim',[-4 4])
    fontsize_figure(gcf,gca,18);
    grid on

%% Load in the SO2 column amount
    %var_UM = 'SO2_perm3_total_column_to_z3000';   %Originally (as of 17th Feb 2021) was using this, but not a
        %good idea since emissions extend up to 4.1km and SO2 gets even
        %higher
    var_UM = 'SO2_perm3_total_column_to_zdomain_top';
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions    
    var_name_out = 'SO2_col_PI_ALL'; 
    clear opts; opts.dummy=NaN; UM_load_var_commands
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out = 'SO2_col_PD_ALL';
    clear opts; opts.dummy=NaN; UM_load_var_commands             
    
    SO2_PI = meanNoNan(SO2_col_PI_ALL,3);
    SO2_PD = meanNoNan(SO2_col_PD_ALL,3);
    
    dSO2_time_mean = meanNoNan(SO2_col_PD_ALL - SO2_col_PI_ALL , 3);
    
    iload_no_orog = 0;
    if iload_no_orog==1
        
        %No orog
        um_case=volc_OFF_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
        var_name_out = 'SO2_col_PI_ALL2';
        clear opts; opts.dummy=NaN; UM_load_var_commands
        
        um_case=volc_ON_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
        var_name_out = 'SO2_col_PD_ALL2';
        clear opts; opts.dummy=NaN; UM_load_var_commands
        
        SO2_PI2 = meanNoNan(SO2_col_PI_ALL2,3);
        SO2_PD2 = meanNoNan(SO2_col_PD_ALL2,3);
        
        dSO2_time_mean2 = meanNoNan(SO2_col_PD_ALL2 - SO2_col_PI_ALL2 , 3);
        
    end
    
     
    iplot_Nd=1;
    if iplot_Nd==1
        %dat_modis = ( SO2_PD - SO2_PI); var_UM = 'Change in SO_2 (PD minus PI; kg m^{-2})';      
        dat_modis = 1e3*dSO2_time_mean; var_UM = 'Change in SO_2 (Volc ON minus OFF)';     
    
        %run plotting script
        figure
        subtitle_str = var_UM;
        ilabel_colorbar = 1; col_bar_lab_str = '(g m^{-2})';
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        clear ilabel_colorbar
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));        
        %caxis([0 300]);
        
        
        dat_modis = dSO2_time_mean2; var_UM = 'No orography change in SO_2 (PD minus PI; kg m^{-2})';     
    
        %run plotting script
        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));        
        %caxis([0 300]);
        
    end
    
    %% Save data in .mat file inc. coarse grained data
    
    save_file = [UM_base_dir um_case_PI '/so2.mat'];
    Hawaii_aod_coarse_grain_save(save_file,SO2_col_PI_ALL,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    
    save_file = [UM_base_dir um_case_PD '/so2.mat'];
    Hawaii_aod_coarse_grain_save(save_file,SO2_col_PD_ALL,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    
%% Load in the sub-grid CF
        
vars_UM = {'CF_lwc_weighted_total_column_to_zdomain_top'};

for ivar=1:length(vars_UM)
    
    var_UM = vars_UM{ivar};
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PI = [var_UM '_PI_ALL'];
    var_name_out = var_name_out_PI;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PD = [var_UM '_PD_ALL'];
    var_name_out = var_name_out_PD;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    eval_str=[var_UM '_PI = meanNoNan(' var_name_out_PI ',3);']; eval(eval_str);
    eval_str=[var_UM '_PD = meanNoNan(' var_name_out_PD ',3);']; eval(eval_str);
    
    eval_str=['d' var_UM ' = meanNoNan(' var_name_out_PD ' - ' var_name_out_PI ',3);']; eval(eval_str);
    
    
    iplot_SW=1;
    if iplot_SW==1
        var_UM = vars_UM{ivar};
        %dat_modis = ( SO2_PD - SO2_PI); var_UM = 'Change in SO_2 (PD minus PI; kg m^{-2})';
        dat_modis = eval(['d' var_UM]); var_UM = ['Change in ' var_UM ' (PD minus PI)'];
        %run plotting script
        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-0.1 0.1]);
        
    end
    
end

%% Load in the sub-grid CF max random overlap
        
vars_UM = {'total_max_random_cloud_amount'};

for ivar=1:length(vars_UM)
    
    var_UM = vars_UM{ivar};
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PI = [var_UM '_PI_ALL'];
    var_name_out = var_name_out_PI;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PD = [var_UM '_PD_ALL'];
    var_name_out = var_name_out_PD;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    eval_str=[var_UM '_PI = meanNoNan(' var_name_out_PI ',3);']; eval(eval_str);
    eval_str=[var_UM '_PD = meanNoNan(' var_name_out_PD ',3);']; eval(eval_str);
    
    eval_str=['d' var_UM ' = meanNoNan(' var_name_out_PD ' - ' var_name_out_PI ',3);']; eval(eval_str);
    
    
    iplot_SW=1;
    if iplot_SW==1
        var_UM = vars_UM{ivar};
        %dat_modis = ( SO2_PD - SO2_PI); var_UM = 'Change in SO_2 (PD minus PI; kg m^{-2})';
        dat_modis = eval(['d' var_UM]); var_UM = ['Change in ' var_UM ' (PD minus PI)'];
        %run plotting script        figure
        subtitle_str = var_UM;
	isave_plot=0;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-0.1 0.1]);
        isave_plot=1; 
	if isave_plot==1
		savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_LWP_' titlenam_driver '_vs_satellite'];
		clear opts
	%	        opts.iplot_png=1;
		opts.iplot_eps=1;
		saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
		%    close(gcf);
	end

        var_UM = vars_UM{ivar};
        icoarse_grain=1;
        M_coarse_grain=10; N_coarse_grain=10;
        %var_UM = '';
        dat_modis = eval(['100*(' var_UM '_PD - '  var_UM '_PI )./ ' var_UM '_PI;']); var_UM = ['% Change in ' var_UM] % \DeltaLWP_{ic subgrid} (g m^{-2})';
        %dat_modis = meanNoNan(W1_orig-W0_orig,3); %similar result
        subtitle_str = var_UM;
        figure
        ilabel_colorbar = 1; col_bar_lab_str = '(%)';
	isave_plot=0;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%        caxis([-50 50]);
        caxis([-20 20]);
        icoarse_grain=0;

	isave_plot=1;
	if isave_plot==1
		savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_LWP_' titlenam_driver '_vs_satellite'];
		clear opts
		        opts.iplot_png=1;
		opts.iplot_eps=1;
		saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
		%    close(gcf);
	end

        
        
        
    end
    
end

isave=1;
if isave==1
    CF_str_var = 'total_max_random_cloud_amount';
    
    save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
    Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PI_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    
    save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
    Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PD_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
end

%% Cloud optical depth estimate
CTT=278; P=850e2; k=0.8;
[tau_PI_ALL,H_PI]=MODIS_tau_func_N_LWP(LWP_PI_ALL*1e-3,Nd_PI_ALL,CTT,k,P); %LWP needs to be in kg/m2 and Nd in per m3
[tau_PD_ALL,H_PD]=MODIS_tau_func_N_LWP(LWP_PD_ALL*1e-3,Nd_PD_ALL,CTT,k,P);


isave=0;
if isave==1
    CF_str_var = 'tau_cloud';
    
    save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
    Hawaii_aod_coarse_grain_save(save_file,eval(['tau_PI_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
    
    save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
    Hawaii_aod_coarse_grain_save(save_file,eval(['tau_PD_ALL']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
end

%% Tau based cloud fraction - better for MODIS comparison?

        thresh_lwp_tau = 0.3;                     
        %thresh_lwp_tau = 1.0; 
        %thresh_lwp_tau = 2.0; 
        cf_str = ['CF_{tau>' num2str(thresh_lwp_tau) '}'];
        CF_str_var = ['CF_tau' num2str(thresh_lwp_tau)]; CF_str_var=remove_character(CF_str_var,'.','pt');
        
        eval_str= [CF_str_var '_PI = zeros(size(LWP_PI_ALL));']; eval(eval_str);
        eval_str= [CF_str_var '_PI(isnan(LWP_PI_ALL)) = NaN;']; eval(eval_str);
        eval_str= [CF_str_var '_PI(tau_PI_ALL>thresh_lwp_tau) = 1;'] ; eval(eval_str);       
        
        eval_str= [CF_str_var '_PD = zeros(size(LWP_PD_ALL));']; eval(eval_str);
        eval_str= [CF_str_var '_PD(isnan(LWP_PD_ALL)) = NaN;']; eval(eval_str);
        eval_str= [CF_str_var '_PD(tau_PD_ALL>thresh_lwp_tau) = 1;'] ; eval(eval_str);          
        
        save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PI']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PD']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
%% LWP based cloud fraction

        thresh_val = 1; %g/m2
        thresh_val = 2; 
        thresh_val = 5; 
        %thresh_val = 10; 
        thresh_val = 20; 
        %thresh_val = 30; 
        
        cf_str = ['CF_{LWP>' num2str(thresh_val) '}'];
        CF_str_var = ['CF_LWP' num2str(thresh_val)]; CF_str_var=remove_character(CF_str_var,'.','pt');
        
        eval_str= [CF_str_var '_PI = zeros(size(LWP_PI_ALL));']; eval(eval_str);
        eval_str= [CF_str_var '_PI(isnan(LWP_PI_ALL)) = NaN;']; eval(eval_str);
        eval_str= [CF_str_var '_PI(LWP_PI_ALL>thresh_val) = 1;'] ; eval(eval_str);       
        
        eval_str= [CF_str_var '_PD = zeros(size(LWP_PD_ALL));']; eval(eval_str);
        eval_str= [CF_str_var '_PD(isnan(LWP_PD_ALL)) = NaN;']; eval(eval_str);
        eval_str= [CF_str_var '_PD(LWP_PD_ALL>thresh_val) = 1;'] ; eval(eval_str);  
        
        
        save_file = [UM_base_dir um_case_PI '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PI']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        save_file = [UM_base_dir um_case_PD '/' CF_str_var '.mat'];
        Hawaii_aod_coarse_grain_save(save_file,eval([CF_str_var '_PD']),gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);
        
        

%% Liquid cloud top height
        
vars_UM = {'Max_Cloud_Height_in_cloud_LWC'};

Hawaii_plot_UM_vars_diffs_etc %Loads the variables and plots PD minus PI diff, PI only and % diff (time
%means).

%% AODs      
vars_UM = {'aod550_aitken','aod550_accum','aod550_coarse','aod550_aitken_insol'}; %Should also have the mineral dust one here, 
%but didn't output that - probably not important.
Hawaii_plot_UM_vars_diffs_etc %Loads the variables and plots PD minus PI diff, PI only and % diff (time
%means).

aod550_total_PI_ALL = aod550_aitken_PI_ALL + aod550_accum_PI_ALL + aod550_coarse_PI_ALL + aod550_aitken_insol_PI_ALL;
aod550_total_PD_ALL = aod550_aitken_PD_ALL + aod550_accum_PD_ALL + aod550_coarse_PD_ALL + aod550_aitken_insol_PD_ALL;

var_UM='aod550_total';
var_UM_clean = remove_character(var_UM ,'.','pt');
var_UM_clean2 = remove_character(var_UM ,'_',' ');
var_name_out_PI = [var_UM_clean '_PI_ALL'];
var_name_out_PD = [var_UM_clean '_PD_ALL'];

eval_str=[var_UM_clean '_PI = meanNoNan(' var_name_out_PI ',3);']; eval(eval_str);
eval_str=[var_UM_clean '_PD = meanNoNan(' var_name_out_PD ',3);']; eval(eval_str);

figure
var_UM = var_UM_clean;
icoarse_grain=1;
M_coarse_grain=10; N_coarse_grain=10;
dat_modis = eval([var_UM '_PD']); var_UM = ['Time mean ' var_UM_clean2 ' (Volcano ON)'];
%run plotting script
subtitle_str = var_UM;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 6]);
icoarse_grain=0;

%% Save AOD data in .mat file inc. coarse grained data
save_file = [UM_base_dir um_case_PI '/aod550_total.mat'];
Hawaii_aod_coarse_grain_save(save_file,aod550_total_PI_ALL,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);

save_file = [UM_base_dir um_case_PD '/aod550_total.mat'];
Hawaii_aod_coarse_grain_save(save_file,aod550_total_PD_ALL,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out_var);

%%
it=14; dat_modis = squeeze(aod550_total_PD_ALL(:,:,it));
figure
var_UM = var_UM_clean;
icoarse_grain=1;
M_coarse_grain=10; N_coarse_grain=10;
var_UM = [var_UM_clean2 ' it=' num2str(it) ' (Volcano ON)'];
%run plotting script
subtitle_str = var_UM;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 6]);
icoarse_grain=0;




%% 3d potemp
vars_UM = {'potential_temp_3d'};


for ivar=1:length(vars_UM)
    
    var_UM = vars_UM{ivar};
    
    var_UM_clean = remove_character(var_UM ,'.','pt');
    var_UM_clean2 = remove_character(var_UM ,'_',' ');
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PI = [var_UM_clean '_PI_ALL'];
    var_name_out = var_name_out_PI;    
%     clear opts; opts.dummy=NaN;
%     opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
%     UM_load_var_commands
    dir_UM = [UM_base_dir um_case_PI];
    dir_UM2 = [dirUM '/' var_UM];
    filename = [dir_UM2 '/merged.nc'];
    nc=netcdf(filename);
    potential_temp_3d_PI_ALL = nc{'potential_temp_3d'}(91:96,:,:,:);
    nc=close(nc);
    
    potential_temp_3d_PI_time_mean = meanNoNan(potential_temp_3d_PI_ALL,1);
    potential_temp_3d_PI_time_mean_prof = meanNoNan(potential_temp_3d_PI_time_mean(:,:),2); 
    
    

    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PD = [var_UM_clean '_PD_ALL'];
    var_name_out = var_name_out_PD;    
%     clear opts; opts.dummy=NaN;
%     opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
%     UM_load_var_commands
    dir_UM = [UM_base_dir um_case_PD];
    dir_UM2 = [dirUM '/' var_UM];
    filename = [dir_UM2 '/merged.nc'];
    nc=netcdf(filename);
    potential_temp_3d_PD_ALL = nc{'potential_temp_3d'}(91:96,:,:,:);
    nc=close(nc);
    
    potential_temp_3d_PD_time_mean = meanNoNan(potential_temp_3d_PD_ALL,1);
    potential_temp_3d_PD_time_mean_prof = meanNoNan(potential_temp_3d_PD_time_mean(:,:),2);
    

end

% Get height info - prob a better way than using the slices...
var_UM = 'SO2_perkg_lon_height_slice_at_ilat=264_lat=19.50';

um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'SO2_slice_PI_ALL';
UM_load_var_commands
z2 = dat_global.gcm_Plat2D_UM(:,1);

figure
hold on
plot(potential_temp_3d_PI_time_mean_prof,z2/1e3,'b');
plot(potential_temp_3d_PD_time_mean_prof,z2/1e3,'r');
legend({'Volc OFF','Volc ON'});
set(gca,'ylim',[0 6]);

%% Mean potemp profiles only for columns that have SO2>thresh
[potential_temp_3d_mean_col_PI,potential_temp_3d_mean_col_PD] = Hawaii_average_profiles_for_high_SO2(UM_base_dir,um_case_PI,um_case_PD,varname_str,SO2_col_PD_ALL)


%% Plot
figure('color','w');
hold on
plot(potential_temp_3d_mean_col_PI,z2/1e3,'b','linewidth',2);
plot(potential_temp_3d_mean_col_PD,z2/1e3,'r','linewidth',2);
legend({'Volc OFF','Volc ON'});
set(gca,'ylim',[0 6]);
xlabel('Potential Temperature (K)');
ylabel('Height (km)');


figure('color','w');
hold on
plot(potential_temp_3d_mean_col_PD - potential_temp_3d_mean_col_PI,z2/1e3,'k','linewidth',2);
legend({'Volc ON minus Volc FFN'});
set(gca,'ylim',[0 6]);
xlabel('Potential Temperature (K)');
ylabel('Height (km)');

%% Mean qV profiles only for columns that have SO2>thresh
varname_str = 'water_vapour_mixing_ratio_3d';
[qV_3d_mean_col_PI,qV_3d_mean_col_PD] = Hawaii_average_profiles_for_high_SO2(UM_base_dir,um_case_PI,um_case_PD,varname_str,SO2_col_PD_ALL);


%% Plot qV profiles
figure('color','w');
hold on
plot(1e3*qV_3d_mean_col_PI,z2/1e3,'b','linewidth',2);
plot(1e3*qV_3d_mean_col_PD,z2/1e3,'r','linewidth',2);
legend({'Volc OFF','Volc ON'});
set(gca,'ylim',[0 6]);
xlabel('Water Vapour Mixing Ratio (g kg^{-1})');
ylabel('Height (km)');


figure('color','w');
hold on
plot(1e3*(qV_3d_mean_col_PD - qV_3d_mean_col_PI),z2/1e3,'k','linewidth',2);
legend({'Volc ON minus Volc FFN'});
set(gca,'ylim',[0 6]);
xlabel('\DeltaWater Vapour Mixing Ratio (g kg^{-1})');
ylabel('Height (km)');
grid on


%% Load in RWP
        
vars_UM = {'RWP_total_column_to_zdomain_top'};

for ivar=1:length(vars_UM)
    
    var_UM = vars_UM{ivar};
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PI = [var_UM '_PI_ALL'];
    var_name_out = var_name_out_PI;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PD = [var_UM '_PD_ALL'];
    var_name_out = var_name_out_PD;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    eval_str=[var_UM '_PI = meanNoNan(' var_name_out_PI ',3);']; eval(eval_str);
    eval_str=[var_UM '_PD = meanNoNan(' var_name_out_PD ',3);']; eval(eval_str);
    
    eval_str=['d' var_UM ' = meanNoNan(' var_name_out_PD ' - ' var_name_out_PI ',3);']; eval(eval_str);
    
    
    iplot_SW=1;
    if iplot_SW==1
        var_UM = vars_UM{ivar};        
        dat_modis = eval(['d' var_UM]); var_UM = ['Change in ' var_UM ' (PD minus PI)'];
        %run plotting script        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-100 100]);
        
        var_UM = vars_UM{ivar};        
        dat_modis = eval(['' var_UM '_PI']); var_UM = [var_UM ' (Volcano OFF)'];
        %run plotting script        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 500]);
        
    end
    
end


%% Load in IWP from sec 30 and ice+snow+graupel from my calculation.
        
vars_UM = {'IWP_sec30','TIWP_total_column_to_zdomain_top'};

for ivar=1:length(vars_UM)
    
    var_UM = vars_UM{ivar};
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PI = [var_UM '_PI_ALL'];
    var_name_out = var_name_out_PI;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PD = [var_UM '_PD_ALL'];
    var_name_out = var_name_out_PD;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    eval_str=[var_UM '_PI = meanNoNan(' var_name_out_PI ',3);']; eval(eval_str);
    eval_str=[var_UM '_PD = meanNoNan(' var_name_out_PD ',3);']; eval(eval_str);
    
    eval_str=['d' var_UM ' = meanNoNan(' var_name_out_PD ' - ' var_name_out_PI ',3);']; eval(eval_str);
    
    
    iplot_SW=1;
    if iplot_SW==1
        var_UM = vars_UM{ivar};        
        dat_modis = eval(['d' var_UM]); var_UM = ['Change in ' var_UM ' (PD minus PI)'];
        %run plotting script        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([-100 100]);
        
        var_UM = vars_UM{ivar};        
        dat_modis = eval(['' var_UM '_PI']); var_UM = [var_UM ' (Volcano OFF)'];
        %run plotting script        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 500]);
        
    end
    
end


%% Load in the SW fluxes
       
%vars_UM = {'SW_up_TOA','SW_up_clear_TOA','SW_up_clean_TOA','SW_up_clean_clear_TOA','SW_down_clean_TOA',...
%    'SW_down_surf','SW_down_clean_surf','SW_down_clean_clear_surf',...
%    'SW_up_clean_surf','SW_up_clean_clear_surf'};

vars_UM = {'SW_up_TOA','SW_up_clean_TOA','SW_up_clean_clear_TOA'};


for ivar=1:length(vars_UM)
    
    var_UM = vars_UM{ivar};
    
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    %um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = 'nested ignore lat lon'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PI = [var_UM '_PI_ALL'];
    var_name_out = var_name_out_PI;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    var_name_out_PD = [var_UM '_PD_ALL'];
    var_name_out = var_name_out_PD;
    clear opts; opts.dummy=NaN;
    opts.time_tol = 3/60/24; %3 min tolerance to allow for the SWTOA timesteps being around 2 mins offset from other vars.
    UM_load_var_commands
    
    eval_str=[var_UM '_PI = meanNoNan(' var_name_out_PI ',3);']; eval(eval_str);
    eval_str=[var_UM '_PD = meanNoNan(' var_name_out_PD ',3);']; eval(eval_str);
    
    eval_str=['d' var_UM ' = meanNoNan(' var_name_out_PD ' - ' var_name_out_PI ',3);']; eval(eval_str);
    
    
    iplot_SW=1;
    if iplot_SW==1
        %dat_modis = ( SO2_PD - SO2_PI); var_UM = 'Change in SO_2 (PD minus PI; kg m^{-2})';
        dat_modis = eval(['d' var_UM]); var_UM = ['Change in ' var_UM ' (PD minus PI; W m^{-2})'];
        
        %run plotting script
        figure
        subtitle_str = var_UM;
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        %caxis([0 300]);
        
    end
    
end

%% Calculate ACI ands ARI forcings
%vars_UM = {'SW_up_TOA','SW_up_clean_TOA','SW_up_clean_clear_TOA','SW_down_clean_TOA',...
%    'SW_down_surf','SW_down_clean_surf','SW_down_clean_clear_surf',...
%    'SW_up_clean_surf','SW_up_clean_clear_surf'};

DRE_ARI_TOA_PI_ALL = - (SW_up_TOA_PI_ALL - SW_up_clean_TOA_PI_ALL);
DRE_ARI_TOA_PD_ALL = - (SW_up_TOA_PD_ALL - SW_up_clean_TOA_PD_ALL); %using -ve since want the forcing in the downwelling direction
ERF_ARI_TOA_ALL = DRE_ARI_TOA_PD_ALL - DRE_ARI_TOA_PI_ALL; 
ERF_ARI_TOA_time_mean = meanNoNan(ERF_ARI_TOA_ALL,3);

CRE_ACI_TOA_PI_ALL = - (SW_up_clean_TOA_PI_ALL - SW_up_clean_clear_TOA_PI_ALL);
CRE_ACI_TOA_PD_ALL = - (SW_up_clean_TOA_PD_ALL - SW_up_clean_clear_TOA_PD_ALL); %using -ve since want the forcing in the downwelling direction
ERF_ACI_TOA_ALL = CRE_ACI_TOA_PD_ALL - CRE_ACI_TOA_PI_ALL; 
ERF_ACI_TOA_time_mean = meanNoNan(ERF_ACI_TOA_ALL,3);

iplot_SW=1;
if iplot_SW==1    
    dat_modis = ERF_ARI_TOA_time_mean; var_UM = ['Time mean ERF_{ARI} (W m^{-2})'];    
    %run plotting script
    figure
    subtitle_str = var_UM;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 300]);
    
    dat_modis = ERF_ACI_TOA_time_mean; var_UM = ['Time mean ERF_{ACI} (W m^{-2})'];    
    %run plotting script
    figure
    subtitle_str = var_UM;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 300]);
    
    dat_modis = meanNoNan(DRE_ARI_TOA_PI_ALL,3); var_UM = ['Time mean DRE_{ARI} Volc OFF (W m^{-2})'];
    %run plotting script
    figure
    subtitle_str = var_UM;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 300]);
    
    dat_modis = meanNoNan(DRE_ARI_TOA_PD_ALL,3); var_UM = ['Time mean DRE_{ARI} Volc ON (W m^{-2})'];
    %run plotting script
    figure
    subtitle_str = var_UM;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 300]);
    
end
    
    %% Timerseries - restrict region to threshold DeltaNd - filtered for each time separately
    thresh_Nd = -1e99; %per cc
    thresh_Nd = 40; %per cc
    thresh_Nd = 10; %per cc
    % NOTE - to get all of the data need to get thresh_Nd to a large
    % negative value (e.g., -1e99);
    dNd = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
    dLWP = LWP_PD_ALL - LWP_PI_ALL;  
    
    title_str = ['\DeltaN_d >' num2str(thresh_Nd) ' cm^{-3}'];
    ylab_str = '\DeltaLWP (g m^{-2})';
    ylims = [-18 105];    
    clear opts; opts.dummy=NaN;
    
    dLWP_timser=Hawaii_plume_timeseries(dLWP,dNd,[thresh_Nd 1e99],time,title_str,ylab_str,ylims,opts);
    
    %% Plot Nd threshold vs dLWP
    dNd = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
    dLWP = LWP_PD_ALL - LWP_PI_ALL; 
        
    thresh_Nd_multi = [-50:10:800];
    %thresh_Nd_multi = [0 50 100];
    clear yvals opts ymean yn ystd
    opts.no_plot=1;
    for i=1:length(thresh_Nd_multi)
        %yvals(:,i)=Hawaii_plume_timeseries(dLWP,dNd,[thresh_Nd_multi(i) 1e99],time,title_str,ylab_str,ylims,opts);        
        imean=find(dNd>=thresh_Nd_multi(i));
        [ymean(i),yn(i),ystd(i)] = meanNoNan(dLWP(imean),1);
    end
    
    %yvals_timemean = meanNoNan(yvals,1);
    figure
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    plot(thresh_Nd_multi,ymean,'linewidth',3);
    xlabel('\DeltaN_d threshold (cm^{-3})');
    ylabel('Time mean \DeltaLWP (g m^{-2})');
    fontsize_figure(gcf,gca,18);
    
    figure
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    plot(thresh_Nd_multi,yn,'linewidth',3);
    %plot(thresh_Nd_multi,ystd./sqrt(yn),'linewidth',3);
    xlabel('\DeltaN_d threshold (cm^{-3})');
    ylabel('N');
    fontsize_figure(gcf,gca,18);
    
%% Plot multiple dLWP for dNd bins with LWP, SO2, longitude and time filters too    
%Using the script to run multiple plots of this.

Hawaii_calc_LWPic_using_subgrid_CF

thresh_LWPs=[-1 5 10 20 30 40 50];
thresh_LWPs=[-1 5 10 20 50];
%thresh_LWPs=[50];
thresh_LWPs=[5];
thresh_LWPs=[-1];

i=1; clear LON_filters
LON_filters{i} = [-160 -155]; i=i+1;
LON_filters{i} = [-165 -160]; i=i+1;
LON_filters{i} = [-170 -165]; i=i+1;
LON_filters{i} = [-175 -170]; i=i+1;
LON_filters{i} = [-175 -160]; i=i+1;
LON_filters{i} = [-361 361]; i=i+1;

i=1; clear plot_types
%plot_types{i} = 'LWP_PD vs Nd_PD'; i=i+1;
%plot_types{i} = 'LWP_PI vs Nd_PD'; i=i+1;
%plot_types{i} = 'dLWP vs Nd_PD'; i=i+1;
%plot_types{i} = 'LWP_PI vs Nd_PI'; i=i+1;
%plot_types{i} = 'dLWP vs dNd'; i=i+1;
plot_types{i} = 'dLWP_subgrid vs dNd'; i=i+1;
plot_types{i} = 'dCF_subgrid vs dNd'; i=i+1;
%plot_types{i} = 'LWP_subgrid PI vs dNd'; i=i+1;
%plot_types{i} = 'LWP_subgrid PD vs dNd'; i=i+1;
%plot_types{i} = 'LWP_subgrid PI vs Nd PI'; i=i+1;
%plot_types{i} = 'LWP_subgrid PD vs Nd PD'; i=i+1;
%plot_types{i} = 'CF_subgrid PD vs Nd PD'; i=i+1;

for iplot_type=1:length(plot_types)
    plot_type = plot_types{iplot_type};
    
    for iLWP=1:length(thresh_LWPs)
        thresh_LWP = thresh_LWPs(iLWP);
        
        for iLON=1:length(LON_filters)
            LON_filter = LON_filters{iLON};
            Hawaii_plots_dLWP_vs_dNd
        end
        
    end
    
end

    
    
           
 %% Plot dLWP for dNd bins with LWP, SO2, longitude and time filters too 
% The idea was to try and restrict the region to the volcanic plume rather
% than including regions where there was no influence at all.  

plot_type = 'dLWP vs dNd';
%plot_type = 'dNd vs dNaccum';

switch plot_type
    case 'dLWP vs dNd'
        var_test = LWP_PD_ALL - LWP_PI_ALL;  
         
        Nd_type = 'CASIM';
        %Nd_type = 'UKCA';
        
        switch Nd_type
            case 'CASIM'
                var_filter = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
            case 'UKCA'
                var_filter = (Nd_UKCA_PD_ALL - Nd_UKCA_PI_ALL)/1e6;
        end
        
        xlab_str = ['\DeltaN_{d ' Nd_type '} (cm^{-3})'];
        y_units_str = '(g m^{-2})';
        ylab_str = ['\DeltaLWP ' y_units_str];
        
    case 'dNd vs dNaccum'
        
        var_test = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
        var_filter = accum_number_z3000_PD_ALL - accum_number_z3000_PI_ALL;
        
        xlab_str = ['\DeltaN_{accum} (m^{-2})'];
        y_units_str = '(cm^{-3})';
        ylab_str = ['\DeltaN_{d ' Nd_type '} ' y_units_str];
        
end
       
%Filtering.
        
    %LWP
    thresh_LWP = 50; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
    thresh_LWP = -1; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
    var_test(inan)=NaN; var_filter(inan)=NaN; %Remove the filter values too to ensure even bin 
        %sampling when using percentiles
        
    thresh_SO2 = -1; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);        
    thresh_SO2 = 1e-5; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);
    %thresh_SO2 = 1e-4; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);
    
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    %CF
    iCF_filter = 0;
    if iCF_filter==1
        thresh_CF = 0.95; inan = find(low_CF_PD_ALL<thresh_CF | low_CF_PI_ALL<thresh_CF);
        thresh_CF = -1; inan = find(low_CF_PD_ALL<thresh_CF | low_CF_PI_ALL<thresh_CF);
        var_test(inan)=NaN; var_filter(inan)=NaN;
        
        CF_filter_str = [', CF threshold = ' num2str(thresh_CF)];
    else
        CF_filter_str = '';
    end
    
    %Longitude
    LON_filter = [-160 -155];
    %LON_filter = [-165 -160];
    %LON_filter = [-170 -165];
    %LON_filter = [-175 -170];
    %LON_filter = [-361 361];    
         
    lon_rep = repmat(gcm_Plon2D_UM,[1 1 size(var_test,3)]);
    inan=find(lon_rep<LON_filter(1) | lon_rep>=LON_filter(2));
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    %Time filtering   
    itime_filter = 0;
    if itime_filter==1
        time_t0 = time_out-time_out(1); %time in days
        time_rep = repmat(time_t0(:),[1 size(var_test,1) size(var_test,2)]);
        time_rep = permute(time_rep,[2 3 1]);
        
        thresh_time = [0 1]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
        %thresh_time = [1 2]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
        %thresh_time = [2 3]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
        thresh_time = [3 4]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
        var_test(inan)=NaN; var_filter(inan)=NaN;
        
        days_str = ['days=' num2str(thresh_time(1)) ' to ' num2str(thresh_time(2))];
    else
        days_str='';
    end
    
    %Make bins using percentiles
    
    %var_filter_bin_edges = [-50:10:800];
    var_filter_bin_edges = prctile(var_filter(:),[0:0.5:100]);
    var_filter_bin_edges = prctile(var_filter(:),[0:1:100]);
        
         
    title_str=[''];
        
    opts=[];                
    [y,std_dev,N]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts); 
    title_str=['LWP threshold = ' num2str(thresh_LWP) ' g m^{-2}, SO_2 threshold = ' num2str(thresh_SO2,'%1.0e')...
        CF_filter_str ', ' num2str(LON_filter(1)) ' to ' num2str(LON_filter(2)) '^{o}E'...
        ', ' days_str...
        ', minN=' num2str(min(N))];       
    
    std_err = std_dev./sqrt(N);    
    yplot = y;
    %yplot(N<=2 | std_err>10 | isnan(std_err)==1)=NaN;    
    %yplot(N<=200 | std_err>10 | isnan(std_err)==1)=NaN;    
    
    figure('color','w');
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
    plot(mid_points,yplot,'bo-','linewidth',3);
    xlabel(xlab_str);
    ylabel(ylab_str);
    title(title_str);
    fontsize_figure(gcf,gca,18);    
    grid on
    
    
    savefile = [savedir_date run_set ' ' title_str];
    %save(savefile,'mid_points','yplot','LON_filter','thresh_LWP','thresh_SO2','thresh_CF','-V7.3');
    
%    loadfile = [savedir_date 'Blending option=3, volcano starting 12 UTC, no orography, adjusted emission height LWP threshold = 50 g m^{-2}, SO_2 threshold = 1e-05, -160 to -155^{o}E, , minN=1372.mat'];
%    dat = load(loadfile);
    
% %% Restrict region to threshold DeltaNd - filtered using time mean
%     thresh_Nd = 100; %per cc
%     dNd_map = (Nd_PD_map - Nd_PI_map)/1e6;
%     %dNd_map2 = repmat(dNd_map
%     
%     
%     iplume = find(dNd_map >= thresh_Nd);  %1d linear indices for 2D map.
%     [ix,iy]=ind2sub(size(dNd_map),iplume);     
%     ix2=repmat(ix,[nT 1]);
%     iy2=repmat(iy,[nT 1]);
%     iz=[1:nT]';
%     iz2=(repmat(iz,[1 length(ix)]))';
%     iz2=iz2(:);   
%     ii = sub2ind(size(dLWP),ix2,iy2,iz2);
%     
%     dLWP_plume = NaN*ones(size(dLWP));    
%     dLWP_plume(ii) = dLWP(ii);    
%     dom_mean_dLWP_plume = meanNoNan(meanNoNan(dLWP_plume,1),1);
%     
%     dNd_map2 = NaN*ones(size(dNd_map));
%     dNd_map2(iplume)=dNd_map(iplume);
%     qpcolor(dNd_map2);
%     caxis([0 400]);
%     title('\DeltaN_d');
%     
%     figure
%     plot(time,dom_mean_dLWP_plume,'linewidth',3);
%     datetick('x','dd');
%   
%     
%     %legend(leg_strs,'location','northwest');
%     xlabel('Time');
%     %ylabel('SW surface forcing (W m^{-2})');
%     ylabel('\DeltaLWP (g m^{-2})');
%     set(gca,'ylim',[-30 50])
%     fontsize_figure(gcf,gca,18);
%     grid on
%     title(['\DeltaN_d thresh=' num2str(thresh_Nd) ' cm^{-3}, filtered with time-av N_d']);    
 

%% Plot dLWP for UKCA dNd bins
    var_filter = (Nd_UKCA_PD_ALL - Nd_UKCA_PI_ALL)/1e6;
    var_test = LWP_PD_ALL - LWP_PI_ALL;  
    %var_filter_bin_edges = [-50:10:800];
    var_filter_bin_edges = prctile(var_filter(:),[0:0.5:100]);
    
    thresh_LWP = 50; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
    var_test(inan)=NaN;        
        
    xlab_str = ['\DeltaN_{d UKCA} (cm^{-3})'];
    y_units_str = '(g m^{-2})';
    ylab_str = ['\DeltaLWP ' y_units_str];
    title_str='';
    opts=[];    
        
    [y,std_dev,N]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts); 
    
    std_err = std_dev./sqrt(N);    
    yplot = y;
    %yplot(yn<=1 | std_err>10 | isnan(std_err)==1)=NaN;    
    yplot(N<=200 | std_err>10 | isnan(std_err)==1)=NaN;    
    
    figure('color','w');
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
    plot(mid_points,yplot,'bo-','linewidth',3);
    xlabel(xlab_str);
    ylabel(ylab_str);
    fontsize_figure(gcf,gca,18);
    title(['LWP threshold = ' num2str(thresh_LWP) ' g m^{-2}']);
    grid on

%% Restrict region to threshold DeltaNd - filtered using time mean - 2    
    thresh_Nd = -1e99; %per cc
    thresh_Nd = 50; %per cc
    dNd_map = (Nd_PD_map - Nd_PI_map)/1e6;
    %dNd_map2 = repmat(dNd_map
    dNd_timeav = repmat(dNd_map,[1 1 size(dNd,3)]);                      

    title_str = ['\DeltaN_d >' num2str(thresh_Nd) ' cm^{-3}, filtered with time-av N_d'];
    ylab_str = '\DeltaLWP (g m^{-2})';
    ylims = [-18 105];    
    clear opts; opts.dummy=NaN;    
    [dLWP_timser_timemean,dLWP_timser_timemean_N,dLWP_timser_timemean_std,inan]=Hawaii_plume_timeseries(dLWP,dNd_timeav,[thresh_Nd 1e99],time,title_str,ylab_str,ylims,opts);
    
    %Do a map of the filtered dNd field used.
    dNd_map2 = dNd_timeav;
    dNd_map2(inan)=NaN;
    qpcolor(dNd_map2(:,:,1));
    caxis([0 400]);
    title('\DeltaN_d');
   
        
%% Plot dLWP for dSO2 bins 

   
    
    var_filter = SO2_col_PD_ALL - SO2_col_PI_ALL;
    var_test = LWP_PD_ALL - LWP_PI_ALL;  
    %Filter out points with low LWP since can still have large dN values
    %and valid LWPs when there are no clouds
    thresh_LWP = 50; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
    var_test(inan)=NaN;
    
    var_filter_bin_edges = [-1:0.2:16]*1e-3;
    var_filter_bin_edges = [-1:0.01:5]*1e-3;
            
    xlab_str = '\DeltaColumn SO_{2 z=3km} (kg m^{-2})';
    y_units_str = '(g m^{-2})';
    ylab_str = ['Time mean \DeltaLWP ' y_units_str];
    title_str='';
    opts=[];
    
        
    [ymean,ystd,yn]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts);
 
    std_err = ystd./sqrt(yn);    
    yplot = ymean;
    %yplot(yn<=1 | std_err>10 | isnan(std_err)==1)=NaN;    
    yplot(yn<=200 | std_err>10 | isnan(std_err)==1)=NaN;    
    
    figure
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
    plot(mid_points,yplot,'linewidth',3);
    xlabel(xlab_str);
    ylabel(ylab_str);
    fontsize_figure(gcf,gca,18);
    grid on
    
    
   
   % 
%% Plot dLWP for dSO2 bins using function + lon bin restriction

    
    thresh_Nd_multi = [-1:0.2:16]*1e-3;
    %thresh_Nd_multi = [0 50 100];
    
    var_filter = SO2_col_PD_ALL - SO2_col_PI_ALL;
    var_test = LWP_PD_ALL - LWP_PI_ALL;                  
    %var_filter_bin_edges = [-1:0.2:16]*1e-3;
    var_filter_bin_edges = [-1:0.01:1]*1e-3;
    
    LON_filter = [-170 -165];
    LAT_filter = [10 18];
    lon_rep = repmat(gcm_Plon2D_UM,[1 1 size(var_test,3)]);
    ilon=find(lon_rep<LON_filter(1) | lon_rep>=LON_filter(2));
    var_test(ilon)=NaN;
            
    xlab_str = '\DeltaColumn SO_{2 z=3km} (kg m^{-2})';
    y_units_str = '(g m^{-2})';
    ylab_str = ['Time mean \DeltaLWP ' y_units_str];
    title_str='';
    opts=[];
    
        
    [ymean,ystd,yn]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts);
 
    std_err = ystd./sqrt(yn);    
    yplot = ymean;
    %yplot(yn<=1 | std_err>10 | isnan(std_err)==1)=NaN;    
    yplot(yn<=200 | std_err>10 | isnan(std_err)==1)=NaN;    
    
    figure
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
    plot(mid_points,yplot,'linewidth',3);
    xlabel(xlab_str);
    ylabel(ylab_str);
    fontsize_figure(gcf,gca,18);
    grid on
    
    
    

%% Calculate in-cloud LWP
%Use total CF or low CF??
thresh_CF_LWPic = 0.05;
LWPic_PI_ALL = LWP_PI_ALL ./ low_CF_PI_ALL;
LWPic_PI_ALL(low_CF_PI_ALL<thresh_CF_LWPic)=NaN;

LWPic_PD_ALL = LWP_PD_ALL ./ low_CF_PD_ALL;
LWPic_PD_ALL(low_CF_PD_ALL<thresh_CF_LWPic)=NaN;
    
%% Load Naccum data

    var_UM = 'accum_number_ukca_total_column_to_z3000';    
    %var_UM2 = 'accum_number_ukca';
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [accum_number_z3000_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    var_UM = 'accum_number_ukca_total_column_to_z3000';    
    %var_UM2 = 'accum_number_ukca';
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [accum_number_z3000_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);         
    
    
%% Plot dLWP for dNaccum bins    
    var_filter = accum_number_z3000_PD_ALL - accum_number_z3000_PI_ALL;
    var_test = LWP_PD_ALL - LWP_PI_ALL;  
        
    thresh_LWP = 50; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
    var_test(inan)=NaN; var_filter(inan)=NaN; %Remove the filter values too to ensure even bin 
        %sampling when using percentiles
    
    thresh = 1e-4; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh);
    %thresh = 1e-5; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh);
    %thresh = 1e-6; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh);
    thresh = -1; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh);
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    LON_filter = [-165 -160];
    LON_filter = [-170 -165];
    LON_filter = [-175 -170];
    LON_filter = [-361 361];
    %LAT_filter = [10 18];
    lon_rep = repmat(gcm_Plon2D_UM,[1 1 size(var_test,3)]);
    inan=find(lon_rep<LON_filter(1) | lon_rep>=LON_filter(2));
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    
    %var_filter_bin_edges = [-50:10:800];
    var_filter_bin_edges = prctile(var_filter(:),[0:0.5:100]);
    var_filter_bin_edges = prctile(var_filter(:),[0:1:100]);
    var_filter_bin_edges = prctile(var_filter(:),[0:2:100]);
    
    opts=[];                
    [y,std_dev,N]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts); 
    
    xlab_str = ['\DeltaN_{accum z3000} (m^{-2})'];
    y_units_str = '(g m^{-2})';
    ylab_str = ['\DeltaLWP ' y_units_str];        
    title_str=['Lon range ' num2str(LON_filter(1)) ' to ' num2str(LON_filter(2)) ' LWP threshold = ' num2str(thresh_LWP) ' g m^{-2}, SO_2 threshold = ' num2str(thresh,'%1.0e') ' minN=' num2str(min(N))];
    
    std_err = std_dev./sqrt(N);    
    yplot = y;
    yplot(N<=2 | std_err>10 | isnan(std_err)==1)=NaN;    
    %yplot(N<=200 | std_err>10 | isnan(std_err)==1)=NaN;    
    
    figure('color','w');
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
    plot(mid_points,yplot,'bo-','linewidth',3);
    xlabel(xlab_str);
    ylabel(ylab_str);
    title(title_str);
    fontsize_figure(gcf,gca,18);    
    grid on
    
    

    %%



    var_thresh = accum_number_z3000_PD_ALL - accum_number_z3000_PI_ALL;
    LWP_PD = LWP_PD_ALL; LWP_PI = LWP_PI_ALL;
    %LWP_PD = LWPic_PD_ALL; LWP_PI = LWPic_PI_ALL;
    
    var_test = LWP_PD - LWP_PI; 
    thresh_LWP = 50; inan = find(LWP_PI<thresh_LWP | LWP_PD<thresh_LWP);
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    thresh = 1e-4; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh);
    thresh = -1; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh);
    thresh = 1e-5; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh);
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    thresh_Nd_multi = [-1:0.04:4]*1e13;
    thresh_Nd_multi = [-1:0.02:8]*1e12;
    thresh_Nd_multi = prctile(var_thresh(:),[0:2:100]);
    
    if min(size(var_thresh) == size(dLWP))==0
       error(['*** DPG - size of var (' num2str(size(dLWP)) ') needs to be the same as the size of var_filter (' num2str(size(var_thresh)) ')***']);
    end
    clear yvals opts ymean yn ystd
    opts.no_plot=1;
    for i=1:length(thresh_Nd_multi)-1
        %yvals(:,i)=Hawaii_plume_timeseries(dLWP,dNd,[thresh_Nd_multi(i) 1e99],time,title_str,ylab_str,ylims,opts);        
        imean=find(var_thresh>=thresh_Nd_multi(i) & var_thresh<thresh_Nd_multi(i+1));
        [ymean(i),yn(i),ystd(i)] = meanNoNan(var_test(imean),1);        
    end
    
    std_err = ystd./sqrt(yn);    
    yplot = ymean;
    %yplot(yn<=200 | std_err>10 | isnan(std_err)==1)=NaN;
    yplot(yn<=2 | std_err>10 | isnan(std_err)==1)=NaN;
    
    %yvals_timemean = meanNoNan(yvals,1);
    figure('color','w');
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( thresh_Nd_multi(1:end-1) + thresh_Nd_multi(2:end) );
    plot(mid_points,yplot,'bo-','linewidth',3);
    xlabel('\DeltaN_{accum z3000} (m^{-2})');
    ylabel('Time mean \DeltaLWP (g m^{-2})');
    title_str=['LWP threshold = ' num2str(thresh_LWP) ' g m^{-2}, SO_2 threshold = ' num2str(thresh,'%1.0e') ' minN=' num2str(min(yn))];
    title(title_str);
    fontsize_figure(gcf,gca,18);
    grid on
    
    figure('color','w');
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    plot(mid_points,yn,'linewidth',3);
    %plot(thresh_Nd_multi,ystd./sqrt(yn),'linewidth',3);
    xlabel('\Deltax ()');
    ylabel('N');
    fontsize_figure(gcf,gca,18);
    grid on
    set(gca,'yscale','log')  
    
    
    figure('color','w');  
    plot(mid_points,std_err,'linewidth',3);
    %plot(thresh_Nd_multi,ystd./sqrt(yn),'linewidth',3);
    xlabel('\Deltax ()');
    ylabel('std dev/\surd(N)');
    fontsize_figure(gcf,gca,18);
    grid on
    set(gca,'yscale','log');
    
    dat_modis = meanNoNan(var_thresh,3);
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = 'Time mean \DeltaN_{accum z3000} (m^{-2})';
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 6e12]);
    icoarse_grain=0;
    
    it=96;
    dat_modis = var_thresh(:,:,it);
    icoarse_grain=0;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['Time mean \DeltaN_{accum z3000} (m^{-2}), it=' num2str(it)];
    var_UM = '';    
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 6e12]);
    icoarse_grain=0;
    
%% Plot dNd for dSO2 bins

    
    thresh_Nd_multi = [-1:0.2:16]*1e-3;
    %thresh_Nd_multi = [0 50 100];
    
    var_thresh = SO2_col_PD_ALL - SO2_col_PI_ALL;      
    var_test = dNd; var_str='N_d (cm^{-3})';
    
    thresh_LWP = 50; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
    var_test(inan)=NaN; var_thresh(inan)=NaN; %Remove the filter values too to ensure even bin 
    
    
    if min(size(var_thresh) == size(var_test))==0
       error(['*** DPG - size of var (' num2str(size(var_test)) ') needs to be the same as the size of var_filter (' num2str(size(var_thresh)) ')***']);
    end
    clear yvals opts ymean yn ystd
    opts.no_plot=1;
    for i=1:length(thresh_Nd_multi)-1
        %yvals(:,i)=Hawaii_plume_timeseries(var_test,dNd,[thresh_Nd_multi(i) 1e99],time,title_str,ylab_str,ylims,opts);        
        imean=find(var_thresh>=thresh_Nd_multi(i) & var_thresh<thresh_Nd_multi(i+1));
        [ymean(i),yn(i),ystd(i)] = meanNoNan(var_test(imean),1);        
    end
    
    std_err = ystd./sqrt(yn);    
    yplot = ymean;
    yplot(yn<=1 | std_err>10 | isnan(std_err)==1)=NaN;
    
    %yvals_timemean = meanNoNan(yvals,1);
    figure
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    mid_points = 0.5*( thresh_Nd_multi(1:end-1) + thresh_Nd_multi(2:end) );
    plot(mid_points,yplot,'linewidth',3);
    xlabel('\DeltaColumn SO_{2 z=3km} (kg m^{-2})');    
    ylabel(['Time mean \Delta' var_str]);
    fontsize_figure(gcf,gca,18);
    grid on
    
    figure
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    plot(mid_points,yn,'linewidth',3);
    %plot(thresh_Nd_multi,ystd./sqrt(yn),'linewidth',3);
    xlabel('\Deltax (cm^{-3})');
    ylabel('N');
    fontsize_figure(gcf,gca,18);
    grid on
    set(gca,'yscale','log')  
    
    
    figure    
    plot(mid_points,ystd./sqrt(yn),'linewidth',3);
    %plot(thresh_Nd_multi,ystd./sqrt(yn),'linewidth',3);
    xlabel('\Deltax (cm^{-3})');
    ylabel('std dev/\surd(N)');
    fontsize_figure(gcf,gca,18);
    grid on
    set(gca,'yscale','log')     
    
    
%% Plot dNd for dNaccum
    var_filter = accum_number_z3000_PD_ALL - accum_number_z3000_PI_ALL;    
    var_test = (Nd_PD_ALL - Nd_PI_ALL)/1e6;   
    var_filter_bin_edges = [-1:0.04:4]*1e13; 
    var_filter_bin_edges = 10.^[12:0.4:13.6];
        
    xlab_str = ['\DeltaN_{accum z3000} (m^{-2})'];
    y_units_str = '(cm^{-3})';
    ylab_str = ['\DeltaN_d CASIM ' y_units_str];
    title_str='';
    opts=[];
    
        
    [y,std_dev,N]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts);
 
    
  %%
    
    
    
    
   
    
    %Plot zonal (over all lons) mean
    [zonal_PI,lats,Nout_PI]=UM_calc_lat_mean(dat_global.gcm_Plat2D_UM,dat_global.gcm_Plon2D_UM,mean_LWP_PI,[],100);
    [zonal_PD,lats,Nout_PD]=UM_calc_lat_mean(dat_global.gcm_Plat2D_UM,dat_global.gcm_Plon2D_UM,mean_LWP_PD,[],100);

    figure
    Nmin=min(Nout_PI,Nout_PD);
    thresh_Nmin = median(Nmin)*0.8;
    i=find(Nmin<thresh_Nmin);
    lats(i)=NaN;
    plot(zonal_PD-zonal_PI,lats);
    
    
    
    
    dat_modis = mean_LWP_PI;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 200]);
    tit_wrapped{3}=um_case_PI;
    title(tit_wrapped);
    
   % Save
    if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'LWP_UKCA'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end    

    dat_modis = mean_LWP_PD;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 200]);
    tit_wrapped{3}=um_case_PD;
    title(tit_wrapped);


   % Save
    if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'LWP_CASIM'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end


%Plot PD minus PI change
    var_UM = var_UM_LWP; 
    dat_modis = mean_LWP_PD - mean_LWP_PI;    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-30 30]);
    tit_wrapped{3}=[um_case_PD ' - ' um_case_PI];
    title(tit_wrapped);

    
    % Save
    if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'LWP_diff'];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end
    

%Plot PD minus PI % change    
    var_UM = 'Prc diff in LWP'
    %dat_modis = 100*2*(mean_LWP_PD - mean_LWP_PI)./(mean_LWP_PI + mean_LWP_PD);  
    dat_modis = 100*(mean_LWP_PD - mean_LWP_PI)./(mean_LWP_PI);      
    dat_modis(mean_LWP_PI<1.0)=NaN;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-100 100]);
    tit_wrapped{3}=['% change ' um_case_PD ' vs ' um_case_PI];
    title(tit_wrapped);

    
    % Save
    if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'LWP_prc_diff'];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end
 
    
    
    
    
    
    
%% Rain rates (LS_surf_rain_rate)
isave_plots=0;

clear gca
 var_UM_LWP = 'LS_surf_rain_rate'; %Might want to try RWP in just BL, rain rate at a particular height, etc.?

    var_UM = var_UM_LWP;
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [LS_rain_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    
    var_UM = var_UM_LWP; 
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [LS_rain_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);   
    
    
    mean_LS_rain_PI = meanNoNan(LS_rain_PI_ALL,3);
    mean_LS_rain_PD = meanNoNan(LS_rain_PD_ALL,3); 
    
 %Plot the mean rain rates 
    dat_modis = mean_LS_rain_PI;
    icoarse_grain=0;
    
    %run plotting script
    figure
    subtitle_str = um_case_PI;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 0.2]);
    %tit_wrapped{3}=um_case_PI;
    %title(tit_wrapped);   
    
       % Save
    if isave_plots==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'LWP_CASIM'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end
    
%Plot the mean rain rates    
    dat_modis = mean_LS_rain_PD;
    icoarse_grain=0;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 0.2]);
    tit_wrapped{3}=um_case_PD;
    title(tit_wrapped);


   % Save
    if isave_plots==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'LWP_CASIM'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end
    
%Plot PD minus PI change

    icoarse_grain=1;
    M_coarse_grain=8; N_coarse_grain=8;

    var_UM = var_UM_LWP; 
    dat_modis = mean_LS_rain_PD - mean_LS_rain_PI;    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-0.03 0.03]);
    tit_wrapped{3}=[um_case_PD ' - ' um_case_PI];
    title(tit_wrapped);

    
    % Save
    if isave_plots==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'LWP_diff'];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end
    

%% Delta LWP combined plot using non-time averaged data loaded from the files
iplot_CASIM_UKCA_combined=1;
if iplot_CASIM_UKCA_combined==1
    %Load both sets of data
    filename_save_dLWP_vs_rain_rates_CASIM = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-bc308/save_dLWP_vs_rain_rates.mat'; 
    filename_save_dLWP_vs_rain_rates_UKCA = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-ba050/save_dLWP_vs_rain_rates.mat'; leg_str{1}='UKCA'; leg_str{2}='CASIM';
    loaded_dat_CASIM=load(filename_save_dLWP_vs_rain_rates_CASIM);
    loaded_dat_UKCA=load(filename_save_dLWP_vs_rain_rates_UKCA);
    
    LS_rain_bins = loaded_dat_CASIM.LS_rain_bins;
    LS_rain_bins_time_mean = loaded_dat_CASIM.LS_rain_bins_time_mean;
    mid_rain_bins = 0.5*(LS_rain_bins(1:end-1) + LS_rain_bins(2:end));
    mid_rain_bins_time_mean = 0.5*(LS_rain_bins_time_mean(1:end-1) + LS_rain_bins_time_mean(2:end));
    
    %Plot mean dLWP vs rain rate bins for non-time av data
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures    
    col_str='b'; h=plot(mid_rain_bins,loaded_dat_UKCA.Y_me,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='UKCA'; set(h,'color',col_str);
    hold on
    col_str='r'; h=plot(mid_rain_bins,loaded_dat_CASIM.Y_me,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('\DeltaLWP (g m^{-2})');
    title('Non time-averaged data');
    legend(leg_str);
    grid on
    
    %Histogram
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures    
    col_str='b'; h=plot(mid_rain_bins,loaded_dat_UKCA.Y_N,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='UKCA'; set(h,'color',col_str);
    hold on
    col_str='r'; h=plot(mid_rain_bins,loaded_dat_CASIM.Y_N,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('No. datapoints');
    title('Non time-averaged data');
    legend(leg_str);
    grid on
    
    %Difference in histograms
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures  
    dN = loaded_dat_CASIM.Y_N - loaded_dat_UKCA.Y_N;
    dat = dN;
    %dat = dN ./ loaded_dat_UKCA.Y_N_time_mean; %normalised by UKCA N values.
    col_str='k'; h=plot(mid_rain_bins,dat,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='CASIM minus UKCA'; set(h,'color',col_str);
    %hold on
    %col_str='r'; h=plot(mid_rain_bins_time_mean,loaded_dat_CASIM.Y_N_time_mean,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    %set(gca,'yscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('Diff in no. datapoints (CASIM minus UKCA)');
    title('Non time-averaged data');
    legend(leg_str);
    grid on

    % Contribution of each bin to the CASIM minus UKCA difference in dLWP
    % where dLWP is the volcano on minus off value
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures   
    N = sum(loaded_dat_UKCA.Y_N); %total number is the same in both CASIM and UKCA cases
    dat = 1/N * ( loaded_dat_CASIM.Y_N.*loaded_dat_CASIM.Y_me - loaded_dat_UKCA.Y_N.*loaded_dat_UKCA.Y_me );    
    col_str='k'; h=plot(mid_rain_bins,dat,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='CASIM minus UKCA'; set(h,'color',col_str);
    %hold on
    %col_str='r'; h=plot(mid_rain_bins_time_mean,loaded_dat_CASIM.Y_N_time_mean,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    %set(gca,'yscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('Contribution (g m^{-2})');
    clear tit_str
    tit_str{1} = 'Contribution of each bin to the';
    tit_str{2} = 'CASIM minus UKCA difference in';
    tit_str{3} = '\DeltaLWP (non time-averaged data)';
    %title('Contribution of each bin to the CASIM minus UKCA difference in \DeltaLWP (time-averaged data)');
    title(tit_str);
    %legend(leg_str);
    grid on
    
%Time averaged data
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures    
    col_str='b'; h=plot(mid_rain_bins_time_mean,loaded_dat_UKCA.Y_me_time_mean,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='UKCA'; set(h,'color',col_str);
    hold on
    col_str='r'; h=plot(mid_rain_bins_time_mean,loaded_dat_CASIM.Y_me_time_mean,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('\DeltaLWP (g m^{-2})');
    title('Time-averaged data');
    legend(leg_str);
    grid on
    
    %Histogram
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures    
    col_str='b'; h=plot(mid_rain_bins_time_mean,loaded_dat_UKCA.Y_N_time_mean,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='UKCA'; set(h,'color',col_str);
    hold on
    col_str='r'; h=plot(mid_rain_bins_time_mean,loaded_dat_CASIM.Y_N_time_mean,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('No. datapoints');
    title('Time-averaged data');
    legend(leg_str);
    grid on
    
    %Difference in histograms
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures  
    dN = loaded_dat_CASIM.Y_N_time_mean - loaded_dat_UKCA.Y_N_time_mean;
    dat = dN;
    %dat = dN ./ loaded_dat_UKCA.Y_N_time_mean; %normalised by UKCA N values.
    col_str='k'; h=plot(mid_rain_bins_time_mean,dat,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='CASIM minus UKCA'; set(h,'color',col_str);
    %hold on
    %col_str='r'; h=plot(mid_rain_bins_time_mean,loaded_dat_CASIM.Y_N_time_mean,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    %set(gca,'yscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('Diff in no. datapoints (CASIM minus UKCA)');
    title('Time-averaged data');
    legend(leg_str);
    grid on
    
    % Contribution of each bin to the CASIM minus UKCA difference in dLWP
    %I.e. the 
    clear leg_str
    figure
    set(gcf,'color','w');
    increase_font_size_map_figures   
    N = sum(loaded_dat_UKCA.Y_N_time_mean); %total number is the same in both CASIM and UKCA cases
    dat = 1/N * ( loaded_dat_CASIM.Y_N_time_mean.*loaded_dat_CASIM.Y_me_time_mean - loaded_dat_UKCA.Y_N_time_mean.*loaded_dat_UKCA.Y_me_time_mean );    
    col_str='k'; h=plot(mid_rain_bins_time_mean,dat,'s-','markerfacecolor',col_str,'linewidth',2); leg_str{1}='CASIM minus UKCA'; set(h,'color',col_str);
    %hold on
    %col_str='r'; h=plot(mid_rain_bins_time_mean,loaded_dat_CASIM.Y_N_time_mean,'o-','markerfacecolor',col_str,'linewidth',2); leg_str{2}='CASIM'; set(h,'color',col_str);
    
    set(gca,'xscale','log');
    %set(gca,'yscale','log');
    xlabel('Rain rate (mm hr^{-1})');
    ylabel('Contribution (g m^{-2})');
    clear tit_str
    tit_str{1} = 'Contribution of each bin to the';
    tit_str{2} = 'CASIM minus UKCA difference in';
    tit_str{3} = '\DeltaLWP (time-averaged data)';
    %title('Contribution of each bin to the CASIM minus UKCA difference in \DeltaLWP (time-averaged data)');
    title(tit_str);
    %legend(leg_str);
    grid on
    
    %N = sum(loaded_dat_UKCA.Y_N_time_mean);
    rain_rate_diff = sum( 1/N*(loaded_dat_CASIM.Y_N_time_mean.*mid_rain_bins_time_mean - loaded_dat_UKCA.Y_N_time_mean.*mid_rain_bins_time_mean));
    %probably better to get it from the actual LS values for volcano OFF
    rain_rate_diff2 = 3600* (meanNoNan(loaded_dat_CASIM.mean_LS_rain_PI(:),1) - meanNoNan(loaded_dat_UKCA.mean_LS_rain_PI(:),1))
    
end
    
    
%% Delta LWP vs surface rain rate (LS_surf_rain_rate)
%copy from Rosenfeld histogram code

iuse_saved_bins=1;
if iuse_saved_bins==1
    filename_loaded_bins = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-bc308/save_dLWP_vs_rain_rates.mat';
    filename_loaded_bins = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-ba050/save_dLWP_vs_rain_rates.mat';
    loaded_bins = load(filename_loaded_bins);
    
    LS_rain_bins = loaded_bins.LS_rain_bins;
    LS_rain_bins_time_mean = loaded_bins.LS_rain_bins_time_mean;
else
    bin_vals = 'choose';
    bin_vals = 'Log10 Percentiles'; nbins=7; deltalog=NaN;
    bin_vals = 'log10'; nbins=NaN; deltalog=0.3;
    
    min_val = 0.0025;
    
    [LS_rain_bins] = UM_Iceland_rain_bins_FUNC(3600*LS_rain_PI_ALL,bin_vals,min_val,nbins,deltalog);
    
end





X_driver = LS_rain_PI_ALL*3600; %convert to mm/hr
X_driver(X_driver<LS_rain_bins(1)) = LS_rain_bins(1)*1.01; %Force values less than a min bin to be = to the min bin so can have a log x scale.
Y_driver = LWP_PD_ALL - LWP_PI_ALL;

Xbins_DRIVER = LS_rain_bins;
Ybins_DRIVER = [-300:50:300];
Ybins_DRIVER = [-100:10:100];

iplot_histo = 0;
if iplot_histo==1
    %Y_driver(inan)=NaN;
    %Z_driver = forcing(inds_PI_clear_low);
    xlabelstr='Surface rain rate (mm hr^{-1})';
    ylabelstr = '\Delta LWP (g m^{-2})';
    
    DRIVER_template_2D_histo_dLWP_vs_rain_rate_Iceland
    shading faceted; %this adds grid lines to the plot - but gives NaNs a
    %colour...
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 100]);
    %label_text_pcolor; %Script to add text labels of the numbers for each block
    title('dLWP vs rain rate histogram');    
end

%Can run the above to get the 2D histogram and then use Y_mean_accurate.
%Or, if don't need the histogram is easier to just do this, which gives the same result:-
clear Y_me Y_N Y_std
for i=1:length(LS_rain_bins)-1
   ime = find(X_driver>=LS_rain_bins(i) & X_driver<LS_rain_bins(i+1));
   %length(ime)
   [Y_me(i),Y_N(i),Y_std(i)] = meanNoNan(Y_driver(ime),1);        
end


mid_rain_bins = 0.5*(LS_rain_bins(1:end-1) + LS_rain_bins(2:end));
figure
plot(mid_rain_bins,Y_me,'bx-');
set(gca,'xscale','log');
xlabel('Rain rate (mm hr^{-1})');
ylabel('\DeltaLWP (g m^{-2})');


%% With time-average values - also saves the data
X_driver_time_mean = meanNoNan(LS_rain_PI_ALL,3)*3600; %convert to mm/hr

if iuse_saved_bins==0    
    bin_vals = 'choose';
    bin_vals = 'Log10 Percentiles'; nbins=7; deltalog=NaN;
    bin_vals = 'log10'; nbins=NaN; deltalog=0.1;
    
    min_val = 0.0025;
    min_val = 0.01;
        
    [LS_rain_bins_time_mean] = UM_Iceland_rain_bins_FUNC(X_driver_time_mean,bin_vals,min_val,nbins,deltalog);
end

X_driver_time_mean(X_driver_time_mean<LS_rain_bins_time_mean(1)) = LS_rain_bins_time_mean(1)*1.01; %Force values less than a min bin to be = to the min bin so can have a log x scale.
Y_driver_time_mean = meanNoNan(LWP_PD_ALL - LWP_PI_ALL , 3);

for i=1:length(LS_rain_bins_time_mean)-1
   ime = find(X_driver_time_mean>=LS_rain_bins_time_mean(i) & X_driver_time_mean<LS_rain_bins_time_mean(i+1));
   [Y_me_time_mean(i),Y_N_time_mean(i),Y_std_time_mean(i)] = meanNoNan(Y_driver_time_mean(ime),1);        
end

mid_rain_bins = 0.5*(LS_rain_bins_time_mean(1:end-1) + LS_rain_bins_time_mean(2:end));
figure
plot(mid_rain_bins,Y_me_time_mean,'bx-');
set(gca,'xscale','log');
xlabel('Rain rate (mm hr^{-1})');
ylabel('\DeltaLWP (g m^{-2})');

dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case_PD];
save_dLWP_vs_rain_rates = [dirUM '/save_dLWP_vs_rain_rates.mat'];

save(save_dLWP_vs_rain_rates,'-V7.3','Y_me_time_mean','Y_N_time_mean','Y_std_time_mean','LS_rain_bins_time_mean');
save(save_dLWP_vs_rain_rates,'-APPEND','Y_me','Y_N','Y_std','LS_rain_bins','mean_LS_rain_PI','mean_LS_rain_PD');





%% RWP, convective LWP, etc.
    
    iload_rwp=0;
    if iload_rwp==1
        var_UM = 'RWP';
        um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
        [RWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);


        var_UM = 'RWP';
        um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
        [RWP_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

        dat_modis = meanNoNan(RWP_PD_ALL,3);

        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        %    caxis([0 1]);

    end

    iload_conv_lwp_rwp=0;
    if iload_conv_lwp_rwp==1
        
        var_UM = 'Conv_LWP';
        um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
   
        %The fisrt time is missing for the convective LWP diags (5-213 etc)
        %for some reason - so try using a special time indices for this, since the time fields are different.

        array_in=[]; %just test to get the indices for now.
        dim=NaN; %don't need dim if just getting the indices
        [out, time_out_conv_LWP, time_inds_conv_LWP, dtime_match_conv_LWP] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);
        
        [Conv_LWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_conv_LWP,load_type,gcm_Plat2D_UM);


        var_UM = 'Conv_LWP';
        um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
        [Conv_LWP_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_conv_LWP,load_type,gcm_Plat2D_UM);

        dat_modis = meanNoNan(Conv_LWP_PD_ALL,3);

        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        caxis([0 150]);
        
        
        var_UM = 'Conv_RWP';
        um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);          
        [Conv_RWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_conv_LWP,load_type,gcm_Plat2D_UM);


        var_UM = 'Conv_RWP';
        um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);
        [Conv_RWP_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_conv_LWP,load_type,gcm_Plat2D_UM);

        dat_modis = meanNoNan(Conv_RWP_PD_ALL,3);

        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        caxis([0 150]);

    end
    
    if iplot_amsre==1
    
    [amsre_tlwp] = amsre_convert_rain_rate(amsre_dat.lwp_amsre, amsre_dat.rain_amsre, amsre_dat.sst_amsre );
    
    lwp_day_night = 'average';
%    lwp_day_night = 'day';
%    lwp_day_night = 'night';
    switch lwp_day_night
        case 'average'
            %average over both day and night
            amsre_LWP_time_period_mean = meanNoNan( meanNoNan(  amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,1:2),4) ,3);
            amsre_TLWP_time_period_mean = meanNoNan( meanNoNan(  amsre_tlwp(:,:,time_inds_AMSRE,1:2),4) ,3);            
        case 'day'
            amsre_LWP_time_period_mean = meanNoNan(  amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,1),3);
            amsre_TLWP_time_period_mean = meanNoNan(  amsre_tlwp(:,:,time_inds_AMSRE,1),3);            
        case 'night'
            amsre_LWP_time_period_mean = meanNoNan(  amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,2),3);
            amsre_TLWP_time_period_mean = meanNoNan(  amsre_tlwp(:,:,time_inds_AMSRE,2),3);            
    end
    
    
    
    %The AMSRE data in this file is upside down and with 0 deg lon at the
    %centre
    [gcm_Plon2D_AMSRE_orig, gcm_Plat2D_AMSRE_orig] = meshgrid([-179.5:179.5],[89.5:-1:-89.5]);
   
    
    %Run subplotting script for UM, CALIPSO and bias
    UM_ACSIS_global_SUBPLOT_commands_LWP
    
    end
    
    
%% Surface rain rate from large-scale    
isave_plot_global_RR=0;

clear gca

 var_UM_RR = 'LS_surf_rain_rate'; %calculated from Python script using QCL after timestep(0-254)
% var_UM = 'LWP_sec30'; %from section 30-405

    var_UM = var_UM_RR;
    um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [RR_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    
    var_UM = var_UM_RR; 
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);   
    [RR_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);   
%Could add an option to check the times here against another run?    
    
    
    mean_RR_PI = meanNoNan(RR_PI_ALL,3);
    mean_RR_PD = meanNoNan(RR_PD_ALL,3);  
    save(save_time_lat_file,'mean_RR_PI','mean_RR_PD','-APPEND','-V7.3');
    
%Calcultate timeseries of mean RR
    dom_mean_RR_PI = meanNoNan(meanNoNan(RR_PI_ALL,1),1);
    dom_mean_RR_PD = meanNoNan(meanNoNan(RR_PD_ALL,1),1);
    
    save(save_time_lat_file,'dom_mean_RR_PI','dom_mean_RR_PD','-APPEND','-V7.3');    
    
    
    dat_modis = mean_RR_PI;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 2e-4]*3600);
    tit_wrapped{3}=um_case_PI;
    title(tit_wrapped);
    
   % Save
    if isave_plot_global_RR==1
        
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/Iceland/' titlenam_driver];
        savename=[savedir_date 'RR_UKCA'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end    
    

    dat_modis = mean_RR_PD;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 2e-4]*3600);
    tit_wrapped{3}=um_case_PD;
    title(tit_wrapped);


   % Save
    if isave_plot_global_RR==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'RR_CASIM'];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end


%Plot PD minus PI change
    dat_modis = mean_RR_PD - mean_RR_PI;    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global   
    caxis([-1e-4 1e-4]*3600);
    tit_wrapped{3}=[um_case_PD ' - ' um_case_PI];
    title(tit_wrapped);

    
    % Save
    if isave_plot_global_RR==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'RR_diff'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end
    
%Plot PD minus PI % change
    var_UM = 'Prc diff in RR'
    dat_modis = 100* (mean_RR_PD - mean_RR_PI) ./ mean_RR_PI;  
    dat_modis(mean_RR_PD<2e-6)=NaN;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global   
    caxis([-20 100]);
    tit_wrapped{3}=['% change ' um_case_PD ' vs ' um_case_PI];
    title(tit_wrapped);

    
    % Save
    if isave_plot_global_RR==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date 'RR_prc_diff'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end


%% Partition the SW_down forcing (PI to PD) into changes due to CF, LWP and
%% SW

%Clear sky transmisison of the atmosphere (SW_down_surf / SW_down_TOA)
transmission_atmos = 1; %model values suggest around 0.8
transmission_atmos = 0.75; %model values suggest around 0.8
transmission_atmos = 0.65; %

% For TOA fluxes also need the clear-sky albedo - can estimate from cloud free regions?
A_surf = 0.15; %Guess for now - value for CF quite sensitive to this (ranging up to 30%), so should prob
  %check what this is from the model
A_surf = 0.12; %value estimated from the 12 noon VOCALS snapshot - N.B. - this uses the estimated transmission since
 %it is the surface albedo assuming no atmosphere above



it_sw = 1:size(SW_down_PI_ALL,3);
    
f0 = low_CF_PI_ALL(:,:,it_sw); f1 = low_CF_PD_ALL(:,:,it_sw);
N0 = Nd_PI_ALL(:,:,it_sw)/1e6; N1 = Nd_PD_ALL(:,:,it_sw)/1e6; %convert to per cc for calc_SW_surf function
W0 = LWP_PI_ALL(:,:,it_sw)./f0; W1 = LWP_PD_ALL(:,:,it_sw)./f1; %Need to supply the in-cloud LWP, so divide by low CF
cf_min = 0.01;
W0(f0<cf_min)=NaN; W1(f1<cf_min)=NaN; %regions with no cloud - set W to NaN since have no cf to divide by
N0(f0<cf_min)=NaN; N1(f1<cf_min)=NaN; %Do the same for Nd for consistency

%Make the values consistent in terms of NaNs since otherwise it can lead to
%inconsistencies when they are averaged. E.g. when calculating SW_estimated_forcing_total_linear
i=find(isnan(N0)==1); W0(i)=NaN; 
i=find(isnan(W0)==1); N0(i)=NaN; 
i=find(isnan(N1)==1); W1(i)=NaN; 
i=find(isnan(W1)==1); N1(i)=NaN; 

%But only NaN the cloud fraction when we are not at zero CF - for zero CF W and N are allowed
%to be NaN.
i=find(isnan(N0)==1 & f0>=cf_min); W0(i)=NaN; f0(i)=NaN;
i=find(isnan(W0)==1 & f0>=cf_min); N0(i)=NaN; f0(i)=NaN;
i=find(isnan(N1)==1 & f1>=cf_min); W1(i)=NaN; f1(i)=NaN;
i=find(isnan(W1)==1 & f1>=cf_min); N1(i)=NaN; f1(i)=NaN;

%Special case where we have no CF in PI, but some cloud in PD - this would
%usually lead to the change in CF being ignored since the PD W and N would
%be NaN. So, set the PI W and N values in theses cases to be equal to the
%PD values, so that this effect is incorporated into the CF change effect
i=find(f0<cf_min & f1>=cf_min); %If they are both <cf_min then they will both be treated as clear-sky in calc_SW_surf
W0_2 = W0; N0_2=N0;
W0_2(i) = W1(i); N0_2(i) = N1(i);

%Do a similar thing for when use PD as the baseline - set W and N to PI
%values if going from zero to some CF between PD and PI
i=find(f1<cf_min & f0>=cf_min);
W1_2 = W1; N1_2=N1;
W1_2(i) = W0(i); N1_2(i) = N0(i);


dat_modis = meanNoNan(W1-W0,3); var_UM = 'Change in in-cloud LWP';
%dat_modis = meanNoNan(100*(N1-N0)./N0,3); var_UM = '% Change in Nd';
%dat_modis = meanNoNan(W1-W0,3) .* meanNoNan(f1-f0,3); var_UM = 'Change in in-cloud LWP x change in CF';
  %Not sure if the above makes sense!
%it=1; dat_modis = meanNoNan(100*(N1(:,:,it)-N0(:,:,it))./N0(:,:,it),3); var_UM = '% Change in Nd';
%it=1; dat_modis = meanNoNan(100*(W1(:,:,it)-W0(:,:,it))./W0(:,:,it),3); var_UM = '% Change in LWPic';
%it=1; dat_modis = meanNoNan(f1(:,:,it)-f0(:,:,it),3); var_UM = 'Change in CF';
%it=1; dat_modis = meanNoNan(SW_PD(:,:,it)-SW_cf2(:,:,it),3); var_UM = 'SW change due to CF';
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-20 10]);


SW_in = SW_in_PD_ALL(:,:,it_sw);
    

%[Ac(i),tau(i),A_f(i),SW_f(i)] = calc_SW_surf(f1,W0,N0,SW_in,A_surf,transmission_atmos);
%[Ac(i),tau(i),A_W(i),SW_W(i)] = calc_SW_surf(f0,W1,N0,SW_in,A_surf,transmission_atmos);
%[Ac(i),tau(i),A_N(i),SW_N(i)] = calc_SW_surf(f0,W0,N1,SW_in,A_surf,transmission_atmos);


switch SW_surf_or_TOA 
    case 'TOA'
            [Ac,tau,T_f,SW_PI] = calc_SW(f0,W0,N0,SW_in,A_surf,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_PD] = calc_SW(f1,W1,N1,SW_in,A_surf,transmission_atmos,cf_min);

            [Ac,tau,T_f,SW_cf] = calc_SW(f1,W0_2,N0_2,SW_in,A_surf,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_lwp] = calc_SW(f0,W1,N0,SW_in,A_surf,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_Nd] = calc_SW(f0,W0,N1,SW_in,A_surf,transmission_atmos,cf_min);

            [Ac,tau,T_f,SW_cf2] = calc_SW(f0,W1_2,N1_2,SW_in,A_surf,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_lwp2] = calc_SW(f1,W0,N1,SW_in,A_surf,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_Nd2] = calc_SW(f1,W1,N0,SW_in,A_surf,transmission_atmos,cf_min);
        
        otherwise
            [Ac,tau,T_f,SW_PI] = calc_SW_surf(f0,W0,N0,SW_in,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_PD] = calc_SW_surf(f1,W1,N1,SW_in,transmission_atmos,cf_min);

            [Ac,tau,T_f,SW_cf] = calc_SW_surf(f1,W0_2,N0_2,SW_in,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_lwp] = calc_SW_surf(f0,W1,N0,SW_in,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_Nd] = calc_SW_surf(f0,W0,N1,SW_in,transmission_atmos,cf_min);

            [Ac,tau,T_f,SW_cf2] = calc_SW_surf(f0,W1_2,N1_2,SW_in,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_lwp2] = calc_SW_surf(f1,W0,N1,SW_in,transmission_atmos,cf_min);
            [Ac,tau,T_f,SW_Nd2] = calc_SW_surf(f1,W1,N0,SW_in,transmission_atmos,cf_min);

end




%% PI to PD difference plots for cf, LWP and Nd

isave_plot_global_LWP=1;
    
% Plot of LWP change PD to PI    - add RWP?
    %Haven't got RWP for PI run yet.
%     
%     LWP_PI_map = meanNoNan(LWP_PI_ALL,3);
%     LWP_PD_map = meanNoNan(LWP_PD_ALL,3);
%     RWP_PI_map = meanNoNan(RWP_PI_ALL,3);
%     RWP_PD_map = meanNoNan(RWP_PD_ALL,3);
    
%     iRWP_diff=0;
%     if iRWP_diff==1
%         dat_modis = (LWP_PD_map + RWP_PD_map)  - (LWP_PI_map + RWP_PI_map); var_UM = 'Change in LWP+RWP (PD minus PI; g m^{-2})';
%     else
% %        dat_modis = LWP_PD_map - LWP_PI_map; var_UM = 'Change in LWP (PD minus PI; g m^{-2})';
%         dat_modis = meanNoNan(W1-W0,3);var_UM = 'Change in in-cloud LWP (PD minus PI; g m^{-2})'; 
%     end
%     %run plotting script
%     figure
%     UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%     lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%     caxis([-10 10]);
%     
%      % Save
%     isave_plot_global_LWP=0;
%     if isave_plot_global_LWP==1
%         savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
%         clear opts
%         %        opts.iplot_png=1;
%         opts.iplot_eps=1;
%         saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
% %        close(gcf);
%     end 
    
    %dat_modis = 100 * (meanNoNan(f1,3)./meanNoNan(f0,3) - 1); var_UM = 'Percentage change in cloud fraction (PD minus PI)';    
    dat_modis = meanNoNan(f1,3) - meanNoNan(f0,3); var_UM = 'Absolute change in cloud fraction (PD minus PI)';     
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-0.1 0.1]);   
    
      % Save
    if isave_plot_global_LWP==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 
    
 dat_modis = 100 * (meanNoNan(W1,3)./meanNoNan(W0,3) - 1); var_UM = 'Percentage change in in-cloud LWP (PD minus PI)';    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-20 20]);   
    
      % Save
    if isave_plot_global_LWP==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 
    
    
    dat_modis = 100 * (meanNoNan(N1,3)./meanNoNan(N0,3) - 1); var_UM = 'Percentage change in Nd (PD minus PI)';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-300 300]);   
    
      % Save

    if isave_plot_global_LWP==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 



%% Do plots of the above paritioning results
isave_plots_partioning=1;
clear gca

isave_plot=0;

switch SW_surf_or_TOA
    case 'TOA';
        var_UM = 'Total estimated TOA SW indirect forcing';
        SW_estimated_forcing_total = meanNoNan(SW_PI,3) - meanNoNan(SW_PD,3);        
    otherwise
        var_UM = 'Total estimated change in surface SW';
        SW_estimated_forcing_total = meanNoNan(SW_PD,3) - meanNoNan(SW_PI,3);
end

dat_modis = meanNoNan(SW_estimated_forcing_total,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-15 15]);


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;

if i_select_region==1
    [iregion_lin,iregion_lin_edges,SW_indirect_estimate_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_indirect_estimate_region_mean = meanNoNan(SW_indirect_estimate_region(:),1);
end

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

var_UM = 'Actual indirect change (forcing) in surface SW';
dat_modis = meanNoNan(indirect_ALL,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-15 15]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
if i_select_region==1
    [iregion_lin,iregion_lin_edges,SW_indirect_actual_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_indirect_actual_region_mean = meanNoNan(SW_indirect_actual_region(:),1)

    SW_estimate_regional_mean_bias = 100* ( 1 - SW_indirect_estimate_region_mean./SW_indirect_actual_region_mean);
end

var_UM = 'Ratio estimated vs actual change';
dat_modis = meanNoNan(SW_estimated_forcing_total,3) ./ meanNoNan(indirect_ALL,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 2]);

var_UM = 'Absolute diff of estimated vs actual change';
dat_modis = meanNoNan(SW_estimated_forcing_total,3) - meanNoNan(indirect_ALL,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([0 2]);

var_UM = 'Percentage difference (estimated vs actual change)';
SW_estimated_forcing = SW_PD-SW_PI;
dat_modis = 100*(meanNoNan(SW_estimated_forcing,3) ./ meanNoNan(indirect_ALL,3) - 1);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

switch SW_surf_or_TOA
    case 'TOA';
        var_UM = 'Change in TOA SW due to CF';
        SW_estimated_forcing_cf = SW_PI - SW_cf;
    otherwise
        var_UM = 'Change in surface SW due to CF';
        SW_estimated_forcing_cf = SW_cf - SW_PI;
end

SW_estimated_forcing_cf_timemean = meanNoNan(SW_estimated_forcing_cf,3);
dat_modis = SW_estimated_forcing_cf_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
if i_select_region==1
    [iregion_lin,iregion_lin_edges,SW_cf_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_cf_region_mean = meanNoNan(SW_cf_region(:),1);
end


switch SW_surf_or_TOA
    case 'TOA'
        var_UM = 'Change in TOA SW forcing due to LWP';
        SW_estimated_forcing_lwp = SW_PI - SW_lwp;
    otherwise        
        var_UM = 'Change in surface SW due to LWP';
        SW_estimated_forcing_lwp = SW_lwp - SW_PI;
        
end
SW_estimated_forcing_lwp_timemean = meanNoNan(SW_estimated_forcing_lwp,3);
dat_modis = SW_estimated_forcing_lwp_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
if i_select_region==1
    [iregion_lin,iregion_lin_edges,SW_lwp_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_lwp_region_mean = meanNoNan(SW_lwp_region(:),1);
end

switch SW_surf_or_TOA
    case 'TOA'
        var_UM = 'Change in TOA SW forcing due to Nd';
        SW_estimated_forcing_Nd = SW_PI - SW_Nd;        
    otherwise
        var_UM = 'Change in surface SW due to Nd';
        SW_estimated_forcing_Nd = SW_Nd - SW_PI;
end

SW_estimated_forcing_Nd_timemean = meanNoNan(SW_estimated_forcing_Nd,3);
dat_modis = SW_estimated_forcing_Nd_timemean; 
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

%restrict to the region of interest
if i_select_region==1
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_Nd_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_Nd_region_mean = meanNoNan(SW_Nd_region(:),1);
end

% ---- IMPORTANT - here it is important to do the meanNoNan averaging first
% since the CF estimate will have a different number of NaNs than the lwp
% and Nd estimates (since for situations with zero CF in the PI, but some
% CF in PD we are using the PD Nd and LWP values).
var_UM = 'Estimated forcing from linear sum';
SW_estimated_forcing_total_linear = meanNoNan(SW_estimated_forcing_cf,3) + meanNoNan(SW_estimated_forcing_lwp,3) + meanNoNan(SW_estimated_forcing_Nd,3);
dat_modis =  SW_estimated_forcing_total_linear;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-15 15]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

if i_select_region==1
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_indirect_estimate_linear_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_indirect_estimate_linear_region_mean = meanNoNan(SW_indirect_estimate_linear_region(:),1)

    SW_estimate_linear_regional_mean_bias = 100* ( 1 - SW_indirect_estimate_linear_region_mean./SW_indirect_actual_region_mean);

end


var_UM = 'Percentage bias estimated linear sum vs actual change';
dat_modis =  100*(1 - meanNoNan(SW_estimated_forcing_total_linear,3) ./ meanNoNan(indirect_ALL,3));
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


var_UM = 'Absolute bias estimated linear sum vs actual change';
dat_modis =  meanNoNan(SW_estimated_forcing_total_linear,3) - meanNoNan(indirect_ALL,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-5 5]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


%I.e., linear sum vs total change, both using SW estimates from the
%calculation.
var_UM = 'Percentage bias estimated linear sum vs estimated total change';
dat_modis =  100*(1 - meanNoNan(SW_estimated_forcing_total_linear,3) ./ meanNoNan(SW_estimated_forcing,3));

%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

if i_select_region==1
    %regional mean
    SW_estimate_regional_linear_mean_bias_vs_total_estimated = 100* ( 1 - SW_indirect_estimate_linear_region_mean./SW_indirect_estimate_region_mean);
end

var_UM = 'Estimated forcing from gridbox mean LWP change';
SW_estimated_forcing_total_linear_LWP = meanNoNan(SW_estimated_forcing_cf,3) + meanNoNan(SW_estimated_forcing_lwp,3);
dat_modis =  SW_estimated_forcing_total_linear_LWP;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-15 15]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


%% Plots showing the proportion of forcing attributable to each variable - PI baseline
SW_tot_linear_timemean_abs = abs(SW_estimated_forcing_cf_timemean) + abs(SW_estimated_forcing_lwp_timemean) + abs(SW_estimated_forcing_Nd_timemean);

var_UM = 'Percentage of abs forcing from \DeltaCF using PI baseline';
dat_modis =  100*abs(SW_estimated_forcing_cf_timemean) ./ SW_tot_linear_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


var_UM = 'Percentage of abs forcing from \DeltaLWP using PI baseline';
dat_modis =  100*abs(SW_estimated_forcing_lwp_timemean) ./ SW_tot_linear_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


var_UM = 'Percentage of abs forcing from \DeltaNd using PI baseline';
dat_modis =  100*abs(SW_estimated_forcing_Nd_timemean) ./ SW_tot_linear_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end



%% ---- using PD with just one for PI switched.

var_UM = 'Change in surface SW due to CF using PD baseline';
SW_estimated_forcing_cf_PD = SW_PD - SW_cf2;
SW_estimated_forcing_cf_PD_timemean = meanNoNan(SW_estimated_forcing_cf_PD,3);
dat_modis = SW_estimated_forcing_cf_PD_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

if i_select_region==1
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_cf_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_cf_region_mean_PD = meanNoNan(SW_cf_region(:),1);
end

var_UM = 'Change in surface SW due to LWP using PD baseline';
SW_estimated_forcing_lwp_PD = SW_PD - SW_lwp2;
SW_estimated_forcing_lwp_PD_timemean = meanNoNan(SW_estimated_forcing_lwp_PD,3);
dat_modis = SW_estimated_forcing_lwp_PD_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

%restrict to the region of interest
if i_select_region==1
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_lwp_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_lwp_region_mean_PD = meanNoNan(SW_lwp_region_PD(:),1);
end


var_UM = 'Change in surface SW due to Nd using PD baseline';
SW_estimated_forcing_Nd_PD = SW_PD - SW_Nd2; 
SW_estimated_forcing_Nd_PD_timemean = meanNoNan(SW_estimated_forcing_Nd_PD,3);
dat_modis = SW_estimated_forcing_Nd_PD_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

%restrict to the region of interest
if i_select_region==1
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_Nd_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_Nd_region_mean_PD = meanNoNan(SW_Nd_region_PD(:),1);
end


% ---- IMPORTANT - here it is important to do the meanNoNan averaging first
% since the CF estimate will have a different number of NaNs than the lwp
% and Nd estimates (since for situations with zero CF in the PD, but some
% CF in PI we are using the PI Nd and LWP values).
var_UM = 'Estimated forcing from linear sum using PD baseline';
SW_estimated_forcing_total_linear_PD = SW_estimated_forcing_cf_PD_timemean + SW_estimated_forcing_lwp_PD_timemean + SW_estimated_forcing_Nd_PD_timemean;
dat_modis =  SW_estimated_forcing_total_linear_PD;

%Note - get different answers depending on how we do the averaging - i.e.,
%average first and then sum, or sum then average - must be due to the NaNs
%in there. E.g. - may not always get a full suite of cf, Nd and W changes.
SW_estimated_forcing_total_linear_PD = meanNoNan(SW_estimated_forcing_cf_PD,3) + meanNoNan(SW_estimated_forcing_lwp_PD,3) + meanNoNan(SW_estimated_forcing_Nd_PD,3);
dat_modis =  SW_estimated_forcing_total_linear_PD;

%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-20 10]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


%restrict to the region of interest
if i_select_region==1
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_indirect_estimate_linear_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_indirect_estimate_linear_region_mean_PD = meanNoNan(SW_indirect_estimate_linear_region_PD(:),1)

    SW_estimate_linear_regional_mean_bias_PD = 100* ( 1 - SW_indirect_estimate_linear_region_mean_PD./SW_indirect_actual_region_mean);

end


var_UM = 'Percentage bias estimated linear sum vs actual change using PD baseline';
dat_modis =  100*(1 - meanNoNan(SW_estimated_forcing_total_linear_PD,3) ./ meanNoNan(indirect_ALL,3));
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

%I.e., linear sum vs total change, both using SW estimates from the
%calculation.
var_UM = 'Percentage bias estimated linear sum vs estimated total change using PD baseline';
dat_modis =  100*(1 - meanNoNan(SW_estimated_forcing_total_linear_PD,3) ./ meanNoNan(SW_estimated_forcing_total,3)); 
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end

%regional mean
SW_estimate_regional_linear_mean_bias_vs_total_estimated_PD = 100* ( 1 - SW_indirect_estimate_linear_region_mean_PD./SW_indirect_estimate_region_mean);


%% Plots showing the proportion of forcing attributable to each variable
SW_tot_linear_PD_timemean_abs = abs(SW_estimated_forcing_cf_PD_timemean) + abs(SW_estimated_forcing_lwp_PD_timemean) + abs(SW_estimated_forcing_Nd_PD_timemean);

var_UM = 'Percentage of abs forcing from \DeltaCF using PD baseline';
dat_modis =  100*abs(SW_estimated_forcing_cf_PD_timemean) ./ SW_tot_linear_PD_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


var_UM = 'Percentage of abs forcing from \DeltaLWP using PD baseline';
dat_modis =  100*abs(SW_estimated_forcing_lwp_PD_timemean) ./ SW_tot_linear_PD_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end


var_UM = 'Percentage of abs forcing from \DeltaNd using PD baseline';
dat_modis =  100*abs(SW_estimated_forcing_Nd_PD_timemean) ./ SW_tot_linear_PD_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

if isave_plots_partioning==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
end





%% List the means
SW_indirect_estimate_region_mean
SW_estimate_regional_mean_bias
SW_cf_region_mean
SW_lwp_region_mean
SW_Nd_region_mean
SW_indirect_estimate_linear_region_mean
SW_estimate_linear_regional_mean_bias
SW_estimate_regional_linear_mean_bias_vs_total_estimated

fprintf(1,'\nUsing PD baseline\n')

SW_cf_region_mean_PD
SW_lwp_region_mean_PD
SW_Nd_region_mean_PD
SW_indirect_estimate_linear_region_mean_PD
SW_estimate_linear_regional_mean_bias_PD
SW_estimate_regional_linear_mean_bias_vs_total_estimated_PD


%% Test plots
it=1;

var_UM = 'Change in surface SW due to CF using PD baseline';
dat_modis = SW_estimated_forcing_cf_PD(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

var_UM = 'Change in surface SW due to LWP using PD baseline';
dat_modis = SW_estimated_forcing_lwp_PD(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

var_UM = 'Change in surface SW due to Nd using PD baseline';
dat_modis = SW_estimated_forcing_Nd_PD(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);


var_UM = 'Surface SW due to CF using PD baseline';
dat_modis = SW_cf2(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-10 10]);

var_UM = 'Surface SW due to LWP using PD baseline';
dat_modis = SW_lwp2(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-10 10]);

var_UM = 'Surface SW due to Nd using PD baseline';
dat_modis = SW_Nd2(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-10 10]);


var_UM = 'Surface SW estimate in PD';
dat_modis = SW_PD(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-10 10]);

%% more test plots
var_UM = 'f';
dat_modis = f0(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));

var_UM = 'W';
dat_modis = W1_2(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));

var_UM = 'Nd';
dat_modis = N1_2(:,:,it);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));

    
    
    
    
    
    
    
%% Prob wont use this for Iceland cases?    
    
    
%% SW forcing vs cloud fraction / height, etc.

% Will just consider clouds in the PI run for now
% consider coarse graining to help make sure met is the same?  
    
    %also look at PD run to look for CF changes
    um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
    
    var_UM = 'low_cloud_amount';
    
%     filename = [dirUM '/' var_UM '/' run_type '_' var_UM '_native res_ALL.mat'];
%     dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);

    
% %    nT = length(eval(['dat_global.' var_UM '_ALL;']));
%     nT = length(time_inds);
%     low_CF_PD_ALL = NaN * ones([size(gcm_Plat2D_UM,1) size(gcm_Plat2D_UM,2) nT]);    
%     it=0;
%     for it_global_diff=1:nT
%         it=it+1;
%         low_CF_PD_ALL(:,:,it) = eval(['dat_global.' var_UM '_ALL{it_global_diff};']);        
%     end
    
%    [low_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);    
%    low_CF_PD_ALL(low_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??    
    
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    
    [iregion_lin,iregion_lin_edges,cf_low]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,low_CF_PI_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_mid]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,mid_CF_PI_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_high]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,high_CF_PI_ALL(:,:,1:end));     
        %CF has one extra time dimension for some reason - is offset by a
        %timestep I think - regrid in time?
        
    [iregion_lin,iregion_lin_edges,cf_low_PD]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,low_CF_PD_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_mid_PD]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,mid_CF_PD_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_high_PD]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,high_CF_PD_ALL(:,:,1:end));                 
            
    [iregion_lin,iregion_lin_edges,cf_low_av]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,low_CF_av_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_mid_av]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,mid_CF_av_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_high_av]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,high_CF_av_ALL(:,:,1:end));             
    
    cf_low_max = max(cf_low,cf_low_PD);
    cf_mid_max = max(cf_mid,cf_mid_PD);
    cf_high_max = max(cf_high,cf_high_PD);    
    
%    forcing_ALL = SW_down_PD_ALL - SW_down_PI_ALL;
    forcing_ALL = indirect_ALL;
    [iregion_lin,iregion_lin_edges,forcing]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,forcing_ALL(:,:,1:end)); 
    
%% Some calculations of forcing in the different low, mid, high cloud combinations 

if icloud_states==1
    
isave_plot_states_bar_chart=0;
isave_plot=0; %for UM_save_plot routine

savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['states_bar'];


    %Using PI clouds
    thresh_CF_states=0.01;
    [forcing_vs_cloud_state_PI,std_forcing_vs_cloud_state_PI,freq_cloud_state_PI,Nii_cloud_state_PI,inds_PI,std_rel_norm_PI,cloud_states_PI] = calc_SW_forcing_cloud_fraction_combinations(forcing,cf_low,cf_mid,cf_high,thresh_CF_states);
    %Using PD clouds
    [forcing_vs_cloud_state_PD,std_forcing_vs_cloud_state_PD,freq_cloud_state_PD,Nii_cloud_state_PD,inds_PD,std_rel_norm_PD,cloud_states_PD] = calc_SW_forcing_cloud_fraction_combinations(forcing,cf_low_PD,cf_mid_PD,cf_high_PD,thresh_CF_states);
    %Using average of PI and PD clouds
    [forcing_vs_cloud_state_AV,std_forcing_vs_cloud_state_AV,freq_cloud_state_AV,Nii_cloud_state_AV,inds_AV,std_rel_norm_AV,cloud_states_AV] = calc_SW_forcing_cloud_fraction_combinations(forcing,cf_low_av,cf_mid_av,cf_high_av,thresh_CF_states);    
    %Using max of PI and PD clouds
    [forcing_vs_cloud_state_max,std_forcing_vs_cloud_state_max,freq_cloud_state_max,Nii_cloud_state_max,inds_max,std_rel_norm_max,cloud_states_max] = calc_SW_forcing_cloud_fraction_combinations(forcing,cf_low_max,cf_mid_max,cf_high_max,thresh_CF_states);        
    %Using PI clouds, but ony when there is no change compared to PD
    cf_thresh_same = 0.000001;
%    cf_thresh_same = 0.01;    
    i_same=find(abs(cf_low-cf_low_PD)<cf_thresh_same & abs(cf_mid-cf_mid_PD)<cf_thresh_same & abs(cf_high-cf_high_PD)<cf_thresh_same );
    Ltot=length(find(isnan(cf_low)==0)); %find total no. points after regional restriction
    frac_left = length(i_same) ./ Ltot
    [forcing_vs_cloud_state_same,std_forcing_vs_cloud_state_same,freq_cloud_state_same,Nii_cloud_state_same,inds_same,std_rel_norm_same,cloud_states_same] = calc_SW_forcing_cloud_fraction_combinations(forcing(i_same),cf_low(i_same),cf_mid(i_same),cf_high(i_same),thresh_CF_states);    
    
    %Creating a combined dataset of the PI and PD clouds with the forcing
    %replicated - hopefully the positive and negative CF chagnes should then
    %cancel out
    forcing_both = cat(3,forcing,forcing);
    cf_low_both = cat(3,cf_low,cf_low_PD);
    cf_mid_both = cat(3,cf_mid,cf_mid_PD);
    cf_high_both = cat(3,cf_high,cf_high_PD);    
    [forcing_vs_cloud_state_both,std_forcing_vs_cloud_state_both,freq_cloud_state_both,Nii_cloud_state_both,inds_both,std_rel_norm_both] = calc_SW_forcing_cloud_fraction_combinations(forcing_both,cf_low_both,cf_mid_both,cf_high_both,thresh_CF_states);
        
    %Overall contribution of each state to total forcing (taking into
    %account frequency of occurence of each state).
    forcing_overall_vs_cloud_state_PI = forcing_vs_cloud_state_PI.*freq_cloud_state_PI ./ sum(freq_cloud_state_PI);
    forcing_overall_vs_cloud_state_PD = forcing_vs_cloud_state_PD.*freq_cloud_state_PD ./ sum(freq_cloud_state_PD);  
    forcing_overall_vs_cloud_state_AV = forcing_vs_cloud_state_AV.*freq_cloud_state_AV ./ sum(freq_cloud_state_AV);  
    forcing_overall_vs_cloud_state_max = forcing_vs_cloud_state_max.*freq_cloud_state_max ./ sum(freq_cloud_state_max);     
    forcing_overall_vs_cloud_state_same = forcing_vs_cloud_state_same.*freq_cloud_state_same ./ sum(freq_cloud_state_same);  
    forcing_overall_vs_cloud_state_both = forcing_vs_cloud_state_both.*freq_cloud_state_both ./ sum(freq_cloud_state_both);     
    
    %Standard error in the overall mean :-
    [me,N,std_dev]=meanNoNan(forcing(:),1);
    std_err = std_dev / sqrt(N);    
    
    frac_cloud_states_PI = freq_cloud_state_PI./sum(freq_cloud_state_PI);
    frac_cloud_states_PD = freq_cloud_state_PD./sum(freq_cloud_state_PD);  
    frac_cloud_states_AV = freq_cloud_state_AV./sum(freq_cloud_state_AV);  
    frac_cloud_states_max = freq_cloud_state_max./sum(freq_cloud_state_max);      
    frac_cloud_states_same = freq_cloud_state_same./sum(freq_cloud_state_same);      
    frac_cloud_states_both = freq_cloud_state_both./sum(freq_cloud_state_both);      
    
    ratio_freq = frac_cloud_states_PD ./ frac_cloud_states_PI;
    %Expressing the change in no. points from PI to PD as a fraction of the
    %total number of (PI) datapoints - so takes into account the
    %frequencies of the different states.
    % Should sum to zero, which nearly does.
    change_as_prc_of_total = 100*(freq_cloud_state_PD - freq_cloud_state_PI)./sum(freq_cloud_state_PI);
    
    %Find the number of indices (gridboxes) that have the same L, M, H CF
    %configuration between PI and PD (although CFs could be different).
    for istate=1:8
        inter(istate)=length(intersect(inds_PD{istate},inds_PI{istate}));
    end
    
    frac_same_state = sum(inter)/sum(freq_cloud_state_PI);
    
    %Fraction of (PI) datapoints for each state that remain in the same
    %state for PD
    frac_same_state_ALL = inter./freq_cloud_state_PI;
    
    dcf_low = cf_low_PD - cf_low;
    dcf_mid = cf_mid_PD - cf_mid;
    dcf_high = cf_high_PD - cf_high;
    
    ilow=find(abs(dcf_low)>0.1); frac_sig_change_low = length(ilow)/Nii_cloud_state_PI;
    imid=find(abs(dcf_mid)>0.1); frac_sig_change_mid = length(imid)/Nii_cloud_state_PI;    
    ihigh=find(abs(dcf_high)>0.1); frac_sig_change_high = length(ihigh)/Nii_cloud_state_PI;  
    
    
    freq_plot_type = 'percentage';
%    freq_plot_type = 'absolute';
    
switch freq_plot_type
    case 'percentage'

        %plot some barcharts of the above
        X=[1:8];
        figure
        bar(X,100*frac_cloud_states_PI);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Frequency of occurrence (%)');
        title('Present day (PI) aerosol run');

        %plot some barcharts of the above
        X=[1:8];
        figure
        bar(X,100*frac_cloud_states_PD);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Frequency of occurrence (%)');
        title('Present day (PD) aerosol run');


        X=[1:8];
        figure
        bar(X,100*frac_cloud_states_AV);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Frequency of occurrence (%)');
        title('PI, PD average CF');
        
        X=[1:8];
        figure
        bar(X,100*frac_cloud_states_both);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Frequency of occurrence (%)');
        title('Both PI and PD');
        

    case 'absolute'

        %plot some barcharts of the above
        X=[1:8];
        figure
        bar(X,freq_cloud_state_PI);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Number of occurrences');
        title('Present day (PI) aerosol run');

        %plot some barcharts of the above
        X=[1:8];
        figure
        bar(X,freq_cloud_state_PD);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Number of occurrences');
        title('Present day (PD) aerosol run');

        X=[1:8];
        figure
        bar(X,freq_cloud_state_AV);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Number of occurrences');
        title('PI, PD average CF');
        
        X=[1:8];
        figure
        bar(X,freq_cloud_state_max);
        increase_font_size_map_figures
        xlabel('Cloud fraction configuration');
        ylabel('Number of occurrences');
        title('Max PI, PD CF');
        
        

end


    

    if isave_plot_states_bar_chart==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/states_bar_chart_freq'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %    close(gcf);
    end
    
        
    X=[1:8];
    figure
    bar(X,100*frac_cloud_states_same);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Frequency of occurrence (%)');
    title('CF for no-change states'); 

    
    if isave_plot_states_bar_chart==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/states_bar_chart_freq'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %    close(gcf);
    end
    
    X=[1:8];
    figure
    bar(X,change_as_prc_of_total);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Percentage change in no. gridpoints PI to PD relative to total in PI (%)');
    title('Frequency change PI to PD');
    
    X=[1:8];
    figure
    bar(X,100*(ratio_freq-1));
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Percentage change in no. gridpoints PI to PD for each state (%)');
    title('Frequency change PI to PD');

    X=[1:8];
    figure
    bar(X,forcing_vs_cloud_state_PI);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Indirect forcing (W m^{-2})');
    title('Mean forcing within each configuration (PI clouds)');
    
    X=[1:8];
    figure
    bar(X,forcing_vs_cloud_state_PD);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Indirect forcing (W m^{-2})');
    title('Mean forcing within each configuration (PD clouds)');
            
    X=[1:8];
    figure
    bar(X,forcing_vs_cloud_state_AV);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Indirect forcing (W m^{-2})');
    title('Mean forcing within each configuration (PI+PD average clouds)');
    
    X=[1:8];
    figure
    bar(X,forcing_vs_cloud_state_max);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Indirect forcing (W m^{-2})');
    title('Mean forcing within each configuration (PI+PD max clouds)');  
    
    X=[1:8];
    figure
    bar(X,forcing_vs_cloud_state_both);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Indirect forcing (W m^{-2})');
    title('Mean forcing within each configuration (both PI+PD clouds)');      
    
    if isave_plot_states_bar_chart==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/states_bar_chart_mean'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %    close(gcf);
    end
    
    X=[1:8];
    figure
    bar(X,forcing_vs_cloud_state_same);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');
    ylabel('Indirect forcing (W m^{-2})');
    title('Mean forcing within each configuration (Clouds with no change)');

    if isave_plot_states_bar_chart==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/states_bar_chart_mean'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %    close(gcf);
    end
    
    ylims=[-3.5 1];
    
    figure
    bar(X,forcing_overall_vs_cloud_state_PI);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');    
    ylabel('Indirect forcing (W m^{-2})');
    title('Conribution to overall forcing (PI clouds)'); 
    set(gca,'ylim',ylims);
    UM_save_plot(gcf,isave_plot,savedir,[save_str '_' get(get(gca,'title'),'string')] );  
    
    
    figure
    bar(X,forcing_overall_vs_cloud_state_PD);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');    
    ylabel('Indirect forcing (W m^{-2})');
    title('Conribution to overall forcing (PD clouds)'); 
    set(gca,'ylim',ylims);
    UM_save_plot(gcf,isave_plot,savedir,[save_str '_' get(get(gca,'title'),'string')] );
    
    figure
    bar(X,forcing_overall_vs_cloud_state_AV);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');    
    ylabel('Indirect forcing (W m^{-2})');
    title('Conribution to overall forcing (PI+PD average clouds)'); 
    set(gca,'ylim',ylims);
    UM_save_plot(gcf,isave_plot,savedir,[save_str '_' get(get(gca,'title'),'string')] );  
    
    figure
    bar(X,forcing_overall_vs_cloud_state_max);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');    
    ylabel('Indirect forcing (W m^{-2})');
    title('Conribution to overall forcing (PI+PD max clouds)'); 
    set(gca,'ylim',ylims);
    UM_save_plot(gcf,isave_plot,savedir,[save_str '_' get(get(gca,'title'),'string')] );
    
    figure
    bar(X,forcing_overall_vs_cloud_state_both);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');    
    ylabel('Indirect forcing (W m^{-2})');
    tit_str = 'Conribution to overall forcing (PI+PD both clouds)'
    title(tit_str);  
    set(gca,'ylim',ylims);    
    UM_save_plot(gcf,isave_plot,savedir,[save_str '_' get(get(gca,'title'),'string')] );
    

    
    
    
    
    
%     if isave_plot_states_bar_chart==1
%         savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/states_bar_chart_tot'];
%         clear opts
%         %        opts.iplot_png=1;
%         opts.iplot_eps=1;
%         saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%         %    close(gcf);
%     end
    
    
   
    figure
    bar(X,forcing_overall_vs_cloud_state_same);
    increase_font_size_map_figures
    xlabel('Cloud fraction configuration');    
    ylabel('Indirect forcing (W m^{-2})');
    title('Conribution to overall forcing (Clouds with no change)');   
    
    if isave_plot_states_bar_chart==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/states_bar_chart_tot'];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %    close(gcf);
    end
    
    
     
    
    % ------------------
end
    
    
    
%% Forcing and histograms by cloud fraction
    %try when looking only at the negative forcings - also removed zero
    %forcings here too
%    forcing(forcing>=0)=NaN;

    tercile_lower = prctile(forcing(:),[100/3]); %lower tercile boundary (prctile ignores NaNs)
%    forcing(forcing>tercile_lower)=NaN;
    tercile_mid = prctile(forcing(:),[200/3]); %lower tercile boundary (prctile ignores NaNs)

    isave_xdat=1;
    isave_xdat_pdf=1;
    clear xdat_save ydat_save xlabelstr_save xdat_save_pdf ydat_save_pdf
    
    %Using PI cloud fraction
    
% Individual low, mid and high plots
    %run 2D histo template script :-    
    X_driver = cf_low;
    Y_driver = forcing;
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save{isave_xdat}=mid_Xbins;
    ydat_save{isave_xdat}=Y_mean;
    xlabelstr_save{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;    
    

    
    %run 2D histo template script :-
    X_driver = cf_mid;
    Y_driver = forcing;
    xlabelstr='Mid cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save{isave_xdat}=mid_Xbins;
    ydat_save{isave_xdat}=Y_mean;
    xlabelstr_save{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;

    axis1D = 'x';
    waterVapourMay2005
    xdat_save_pdf{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;


    
    %run 2D histo template script :-
    X_driver = cf_high;
    Y_driver = forcing;
    xlabelstr='High cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save{isave_xdat}=mid_Xbins;
    ydat_save{isave_xdat}=Y_mean;
    xlabelstr_save{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;
    
    axis1D = 'x';
    waterVapourMay2005
    xdat_save_pdf{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;
    
%% Plot the mean SW vs cloud fraction bin (low, mid and high on same plot)
isave_plot_global_mean_SW=0;

lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save{iplot},ydat_save{iplot},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{iplot},' fraction','');
hold on

iplot=iplot+1;
plot(xdat_save{2},ydat_save{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{iplot},' fraction','');

iplot=iplot+1;
plot(xdat_save{3},ydat_save{3},'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{iplot},' fraction','');

legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
%ylabel('SW surface forcing (W m^{-2})');
ylabel('SW surface indirect forcing (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);

   


if isave_plot_global_mean_SW==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/SW_vs_CF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    

%% Plot the PDF of cloud fraction for low, mid and high clouds
isave_plot_global_pdf=0;

lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','');
hold on

iplot=iplot+1;
plot(xdat_save_pdf{2},ydat_save_pdf{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','')

iplot=iplot+1;
plot(xdat_save_pdf{3},ydat_save_pdf{3},'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','')

legend(leg_strs,'location','northeast');
xlabel('Cloud fraction');
ylabel('No. data points');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  

if isave_plot_global_pdf==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/CF_PDF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    


%Back up the saved arrays and label as PI arrays
xdat_save_PI = xdat_save;
ydat_save_PI = ydat_save;
xlabelstr_save_PI = xlabelstr_save;

xdat_save_pdf_PI = xdat_save_pdf;
ydat_save_pdf_PI = ydat_save_pdf;
xlabelstr_save_pdf_PI = xlabelstr_save_pdf;


%% Plot the constribution of each cloud fraction to total SW forcing for
%% low, mid and high clouds
isave_plot_global_cont=0;

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf{iplot});
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot}.*ydat_save{iplot}/N,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','');
hold on

iplot=iplot+1;
N = sum(ydat_save_pdf{iplot});
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot}.*ydat_save{iplot}/N,'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','')

iplot=iplot+1;
N = sum(ydat_save_pdf{iplot});
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot}.*ydat_save{iplot}/N,'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','');

legend(leg_strs,'location','southeast');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  


if isave_plot_global_cont==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end 


%% Same, but using PD cloud fractions
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    %run 2D histo template script :-    
    X_driver = cf_low_PD;
    Y_driver = forcing;
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save{isave_xdat}=mid_Xbins;
    ydat_save{isave_xdat}=Y_mean;
    xlabelstr_save{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;    
    

    
    %run 2D histo template script :-
    X_driver = cf_mid_PD;
    Y_driver = forcing;
    xlabelstr='Mid cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save{isave_xdat}=mid_Xbins;
    ydat_save{isave_xdat}=Y_mean;
    xlabelstr_save{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;

    axis1D = 'x';
    waterVapourMay2005
    xdat_save_pdf{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;


    
    %run 2D histo template script :-
    X_driver = cf_high_PD;
    Y_driver = forcing;
    xlabelstr='High cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save{isave_xdat}=mid_Xbins;
    ydat_save{isave_xdat}=Y_mean;
    xlabelstr_save{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;
    
    axis1D = 'x';
    waterVapourMay2005
    xdat_save_pdf{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;
    
    
%% Repeat, but using PI+PD average cloud fractions
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    %run 2D histo template script :-    
    X_driver = cf_low_av;
    Y_driver = forcing;
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av{isave_xdat}=mid_Xbins;
    ydat_save_av{isave_xdat}=Y_mean;
    xlabelstr_save_av{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_av{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;    
    

    
    %run 2D histo template script :-
    X_driver = cf_mid_av;
    Y_driver = forcing;
    xlabelstr='Mid cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av{isave_xdat}=mid_Xbins;
    ydat_save_av{isave_xdat}=Y_mean;
    xlabelstr_save_av{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;

    axis1D = 'x';
    waterVapourMay2005
    xdat_save_pdf_av{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;


    
    %run 2D histo template script :-
    X_driver = cf_high_av;
    Y_driver = forcing;
    xlabelstr='High cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av{isave_xdat}=mid_Xbins;
    ydat_save_av{isave_xdat}=Y_mean;
    xlabelstr_save_av{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;
    
    axis1D = 'x';
    waterVapourMay2005
    xdat_save_pdf_av{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;
    
%% Repeat again, using PI+PD average cloud fractions, but restricting to
%% the low only cloud state
    clear gca
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    %run 2D histo template script :-    
    X_driver = cf_low_av(inds_AV{2});
    Y_driver = forcing(inds_AV{2});
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av_low_only{isave_xdat}=mid_Xbins;
    ydat_save_av_low_only{isave_xdat}=Y_mean;
    xlabelstr_save_av_low_only{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_av_low_only{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av_low_only{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av_low_only{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf_low_only=isave_xdat_pdf+1;    
    
%% Repeat again, using PI cloud fractions, but restricting to
%% the low only cloud state
    clear gca
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    %run 2D histo template script :-    
    X_driver = cf_low(inds_PI{2});
    Y_driver = forcing(inds_PI{2});
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av_low_only_PI{isave_xdat}=mid_Xbins;
    ydat_save_av_low_only_PI{isave_xdat}=Y_mean;
    xlabelstr_save_av_low_only_PI{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_av_low_only_PI{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av_low_only_PI{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av_low_only_PI{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf_low_only=isave_xdat_pdf+1;    
    
    
%% Repeat again, using PD cloud fractions, but restricting to
%% the low only cloud state
    clear gca
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    %run 2D histo template script :-    
    X_driver = cf_low_PD(inds_PD{2});
    Y_driver = forcing(inds_PD{2});
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av_low_only_PD{isave_xdat}=mid_Xbins;
    ydat_save_av_low_only_PD{isave_xdat}=Y_mean;
    xlabelstr_save_av_low_only_PD{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_av_low_only_PD{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av_low_only_PD{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av_low_only_PD{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf_low_only=isave_xdat_pdf+1;    

%% Repeat again, using PI cloud fractions, but restricting to
%% the low only cloud state combined with the clear state
    clear gca
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    %run 2D histo template script :-    
    X_driver = cf_low(inds_PI{2});
    Y_driver = forcing(inds_PI{2});
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av_low_only_PI{isave_xdat}=mid_Xbins;
    ydat_save_av_low_only_PI{isave_xdat}=Y_mean;
    xlabelstr_save_av_low_only_PI{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_av_low_only_PI{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av_low_only_PI{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av_low_only_PI{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1; 
    

%% Repeat again, using PI cloud fractions, but restricting to
%% the low only cloud state combined with the clear state
    clear gca
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    inds_PI_clear_low = [inds_PI{1}; inds_PI{2}]; %Both states 1 and 2 (clear and low-only)
    
    %run 2D histo template script :-    
    X_driver = cf_low(inds_PI_clear_low);
    Y_driver = forcing(inds_PI_clear_low);
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_clearlow_only_PI{isave_xdat}=mid_Xbins;
    ydat_save_clearlow_only_PI{isave_xdat}=Y_mean;
    xlabelstr_save_clearlow_only_PI{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_clearlow_only_PI{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_clearlow_only_PI{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_clearlow_only_PI{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;    
    
    
    
%% Repeat again, using PD cloud fractions, but restricting to
%% the low only cloud state combined with the clear state
    clear gca
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    inds_PD_clear_low = [inds_PD{1}; inds_PD{2}]; %Both states 1 and 2 (clear and low-only)
    
    %run 2D histo template script :-    
    X_driver = cf_low_PD(inds_PD_clear_low);
    Y_driver = forcing(inds_PD_clear_low);
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_clearlow_only_PD{isave_xdat}=mid_Xbins;
    ydat_save_clearlow_only_PD{isave_xdat}=Y_mean;
    xlabelstr_save_clearlow_only_PD{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_clearlow_only_PD{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_clearlow_only_PD{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_clearlow_only_PD{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf=isave_xdat_pdf+1;    
    
    
    
    
%% Repeat again, using low cloud fractions only when there is no change between PI and PD, but restricting to
%% the low only cloud state
    clear gca
    isave_xdat = 1;
    isave_xdat_pdf = 1;
    
    
    %run 2D histo template script :-    
    X_driver = cf_low(i_same(inds_same{2})); 
    Y_driver = forcing(i_same(inds_same{2}));
    xlabelstr='Low cloud fraction';
    DRIVER_template_2D_PDF_02_SW_forcing_vs_CF_UM    
    plot(mid_Xbins,Y_mean,'wo');
    set(gca,'ylim',[-90 10]);
    xdat_save_av_low_only_same{isave_xdat}=mid_Xbins;
    ydat_save_av_low_only_same{isave_xdat}=Y_mean;
    xlabelstr_save_av_low_only_same{isave_xdat} = xlabelstr;
    isave_xdat=isave_xdat+1;    
    
    axis1D = 'x';
    waterVapourMay2005;
    xdat_save_pdf_av_low_only_same{isave_xdat_pdf}=xdat(1).x;
    ydat_save_pdf_av_low_only_same{isave_xdat_pdf}=ydat(1).y;
    xlabelstr_save_pdf_av_low_only_same{isave_xdat_pdf} = xlabelstr;
    isave_xdat_pdf_low_only=isave_xdat_pdf+1;    
    

    
%% 2D mean contribution to forcing from different combinations of PI and PD
%% cloud fraction using PI cloud fractions - start with no resctrictions to
%% states? Or use combined states from PI and PD?
    clear gca
    isave_mean_2D_PDF=1;
    
    
    %inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    inds_PI_clear_low = unique( [inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    
    %run 2D histo template script :-    
    X_driver = cf_low(inds_PI_clear_low);
    Y_driver = cf_low_PD(inds_PI_clear_low);    
    Z_driver = forcing(inds_PI_clear_low);
    xlabelstr='Low cloud fraction Pre-Industrial';
    ylabelstr = 'Low cloud fraction Present Day';

    DRIVER_template_2D_mean_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]);
    title('Mean forcing (W m^{-2})');
    save_str = 'forcing_mean_2D_PDF_cf_PI_PD';
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  
        
    DRIVER_template_2D_overall_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis([-2 2]);    
    %caxis([-2.5 2.5]);       
    title('Overall forcing (W m^{-2})');
    save_str = 'forcing_overall_2D_PDF_cf_PI_PD';
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  

    
    %2D plots of forcing vs cloud state PI and PD
    X_driver = cloud_states_PI;
    Y_driver = cloud_states_PD;    
    Z_driver = forcing;
    xlabelstr='Cloud State Pre-Industrial';
    ylabelstr = 'Cloud State Present Day';
    
    DRIVER_template_2D_overall_SW_forcing_vs_cloudstate_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis([-2 2]);    
    %caxis([-2.5 2.5]);       
    title('Overall forcing (W m^{-2})');
    set(gca,'xlim',[0.5 8.5]);
    set(gca,'ylim',[0.5 8.5]);    
    save_str = 'forcing_overall_2D_PDF_cf_PI_PD';
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  
    
%% Plot the mean SW vs cloud fraction bin
isave_plot_global_cont=0;

lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save{iplot},ydat_save{iplot},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{iplot},' fraction','');
hold on

iplot=iplot+1;
plot(xdat_save{2},ydat_save{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{iplot},' fraction','');

iplot=iplot+1;
plot(xdat_save{3},ydat_save{3},'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{iplot},' fraction','');

legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
%ylabel('SW surface forcing (W m^{-2})');
ylabel('SW surface indirect forcing (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);

   


if isave_plot_global_mean_SW==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/SW_vs_CF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    



%% Plot again, but using average of PI and PD cloud fractions.
isave_plot_global_cont=0;

lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save_av{1},ydat_save_av{1},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_av{1},' fraction',' AV');
hold on

iplot=iplot+1;
plot(xdat_save_av{2},ydat_save_av{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_av{2},' fraction',' AV');

iplot=iplot+1;
plot(xdat_save_av{3},ydat_save_av{3},'k','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_av{3},' fraction',' AV');

legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
%ylabel('SW surface forcing (W m^{-2})');
ylabel('SW surface indirect forcing (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);

   


if isave_plot_global_mean_SW==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/SW_vs_CF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    

%% Plot again using average of PI and PD cloud fractions, but for low CF
%% state only
isave_plot_global_cont=0;

lwidth=4;
fsize=18;
clear leg_strs

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save_av_low_only{1},ydat_save_av_low_only{1},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_av_low_only{1},' fraction',' AV');
hold on


legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
%ylabel('SW surface forcing (W m^{-2})');
ylabel('SW surface indirect forcing (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);

   


if isave_plot_global_mean_SW==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/SW_vs_CF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    




%% Plot with PI and PD on the same plot
isave_plot_global_mean_SW=0;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save_PI{iplot},ydat_save_PI{iplot},'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_PI{iplot},' fraction',' PI');
hold on

iplot=iplot+1;
plot(xdat_save_PI{2},ydat_save_PI{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_PI{iplot},' fraction',' PI');

iplot=iplot+1;
plot(xdat_save_PI{3},ydat_save_PI{3},'k','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_PI{iplot},' fraction',' PI');

%PD clouds
iplot=iplot+1;
plot(xdat_save{1},ydat_save{1},'b--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{1},' fraction',' PD');
hold on

iplot=iplot+1;
plot(xdat_save{2},ydat_save{2},'r--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{2},' fraction',' PD');

iplot=iplot+1;
plot(xdat_save{3},ydat_save{3},'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save{3},' fraction',' PD');

legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
%ylabel('SW surface forcing (W m^{-2})');
ylabel('SW surface indirect forcing (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);

   


if isave_plot_global_mean_SW==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/SW_vs_CF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    





%% Plot the PDF of cloud fraction for just PD clouds
isave_plot_global_pdf=0;

lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','');
hold on

iplot=iplot+1;
plot(xdat_save_pdf{2},ydat_save_pdf{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','')

iplot=iplot+1;
plot(xdat_save_pdf{3},ydat_save_pdf{3},'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','')

legend(leg_strs,'location','northeast');
xlabel('Cloud fraction');
ylabel('No. data points');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  

if isave_plot_global_pdf==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/CF_PDF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    




%% Plot the PDF of cloud fraction for average of PI and PD clouds
isave_plot_global_pdf=0;

lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save_pdf_av{iplot},ydat_save_pdf_av{iplot},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av{iplot},' fraction','');
hold on

iplot=iplot+1;
plot(xdat_save_pdf_av{2},ydat_save_pdf_av{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av{iplot},' fraction','')

iplot=iplot+1;
plot(xdat_save_pdf_av{3},ydat_save_pdf_av{3},'k','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av{iplot},' fraction','')

legend(leg_strs,'location','northeast');
xlabel('Cloud fraction');
ylabel('No. data points');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  

if isave_plot_global_pdf==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/CF_PDF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    

%% Plot the PDF of cloud fraction for average of PI and PD clouds for only low-only cloud state
isave_plot_global_pdf=0;

lwidth=4;
fsize=18;

figure
clear leg_strs
iplot=0;

iplot=iplot+1;
plot(xdat_save_pdf_av_low_only{iplot},ydat_save_pdf_av_low_only{iplot},'linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av_low_only{iplot},' fraction','');
hold on


legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
ylabel('No. data points');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  

if isave_plot_global_pdf==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/CF_PDF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end    


%% Plot the PDF of cloud fraction for PI and PD clouds separately
isave_plot_global_pdf=0;

lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
plot(xdat_save_pdf_PI{iplot},ydat_save_pdf_PI{iplot},'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_PI{iplot},' fraction',' PI');
hold on

iplot=iplot+1;
plot(xdat_save_pdf_PI{2},ydat_save_pdf_PI{2},'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_PI{iplot},' fraction',' PI')

iplot=iplot+1;
plot(xdat_save_pdf_PI{3},ydat_save_pdf_PI{3},'k','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_PI{iplot},' fraction',' PI')

%PD
iplot=iplot+1;
plot(xdat_save_pdf{1},ydat_save_pdf{1},'b--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{1},' fraction',' PD');
hold on

iplot=iplot+1;
plot(xdat_save_pdf{2},ydat_save_pdf{2},'r--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{2},' fraction',' PD')

iplot=iplot+1;
plot(xdat_save_pdf{3},ydat_save_pdf{3},'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{3},' fraction',' PD')

legend(leg_strs,'location','northeast');
xlabel('Cloud fraction');
ylabel('No. data points');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  

if isave_plot_global_pdf==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/CF_PDF_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end 


%% Plot the contribution of each cloud fraction to total SW forcing just for PD
isave_plot_global_cont=0

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf{iplot});
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot}.*ydat_save{iplot}/N,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','');
hold on

iplot=iplot+1;
N = sum(ydat_save_pdf{iplot});
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot}.*ydat_save{iplot}/N,'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','')

iplot=iplot+1;
N = sum(ydat_save_pdf{iplot});
plot(xdat_save_pdf{iplot},ydat_save_pdf{iplot}.*ydat_save{iplot}/N,'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{iplot},' fraction','');

legend(leg_strs,'location','southeast');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  


if isave_plot_global_cont==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end 

%% Plot the contribution of each cloud fraction to total SW forcing using
%% PI and PD cloud fraction average
isave_plot_global_cont=0;

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf_av{iplot});
plot(xdat_save_pdf_av{iplot},ydat_save_pdf_av{iplot}.*ydat_save_av{iplot}/N,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av{iplot},' fraction',' AV');
hold on

iplot=iplot+1;
N = sum(ydat_save_pdf_av{iplot});
plot(xdat_save_pdf_av{iplot},ydat_save_pdf_av{iplot}.*ydat_save_av{iplot}/N,'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av{iplot},' fraction',' AV')

iplot=iplot+1;
N = sum(ydat_save_pdf_av{iplot});
plot(xdat_save_pdf_av{iplot},ydat_save_pdf_av{iplot}.*ydat_save_av{iplot}/N,'k','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av{iplot},' fraction',' AV');

legend(leg_strs,'location','southeast');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  


if isave_plot_global_cont==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end 


%% Plot the contribution of each cloud fraction to total SW forcing using
%% PI and PD cloud fraction average - but only for low-only cloud state
isave_plot_global_cont=0;

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf_av_low_only{iplot});
plot(xdat_save_pdf_av_low_only{iplot},ydat_save_pdf_av_low_only{iplot}.*ydat_save_av_low_only{iplot}/N,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av_low_only{iplot},' fraction',' AV');
hold on


legend(leg_strs,'location','southwest');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  


if isave_plot_global_cont==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end 



%% Plot the contribution of each cloud fraction to total SW forcing using
%% PD cloud fraction average - but only for low-only cloud state
isave_plot_global_cont=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['Forcing_vs_CF'];

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf_av_low_only{iplot});
y_clear_add = zeros(size(ydat_save_pdf_av_low_only{iplot})); y_clear_add(1) = forcing_overall_vs_cloud_state_PD(1);

%y_dat_SW_overall_cont_PD = ydat_save_pdf_av_low_only_PD{iplot}.*ydat_save_av_low_only_PD{iplot}/N;
y_dat_SW_overall_cont_PD = y_clear_add + ydat_save_pdf_av_low_only_PD{iplot}.*ydat_save_av_low_only_PD{iplot}/N;
plot(xdat_save_pdf_av_low_only{iplot},y_dat_SW_overall_cont_PD,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av_low_only{iplot},' fraction',' PD');
hold on


legend(leg_strs,'location','southwest');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
UM_save_plot(gcf,isave_plot_global_cont,savedir,[save_str '_' leg_strs{iplot}] );  


% if isave_plot_global_cont==1
%     savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
%     clear opts
%     %        opts.iplot_png=1;
%     opts.iplot_eps=1;
%     saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
% %    close(gcf);
% end 

%% Plot the contribution of each cloud fraction to total SW forcing using
%% PI cloud fractions - but only for low-only cloud state
isave_plot_global_cont=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['Forcing_vs_CF'];

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf_av_low_only{iplot});
y_clear_add = zeros(size(ydat_save_pdf_av_low_only{iplot})); y_clear_add(1) = forcing_overall_vs_cloud_state_PI(1);
%y_dat_SW_overall_cont_PI = ydat_save_pdf_av_low_only_PI{iplot}.*ydat_save_av_low_only_PI{iplot}/N;
y_dat_SW_overall_cont_PI = y_clear_add + ydat_save_pdf_av_low_only_PI{iplot}.*ydat_save_av_low_only_PI{iplot}/N;
plot(xdat_save_pdf_av_low_only{iplot},y_dat_SW_overall_cont_PI,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av_low_only{iplot},' fraction',' PI');
hold on


legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
UM_save_plot(gcf,isave_plot_global_cont,savedir,[save_str '_' leg_strs{iplot}] );  


% if isave_plot_global_cont==1
%     savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
%     clear opts
%     %        opts.iplot_png=1;
%     opts.iplot_eps=1;
%     saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
% %    close(gcf);
% end 

%% PI and PD for the above on the same plot
isave_plot_global_cont=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['Forcing_vs_CF'];

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0; clear leg_strs

iplot=iplot+1;
plot(xdat_save_pdf_av_low_only{1},y_dat_SW_overall_cont_PI,'b','linewidth',lwidth);
%leg_strs{iplot} = remove_character(xlabelstr_save_pdf_av_low_only{1},' fraction',' PI');
leg_strs{iplot} = 'Pre-Industrial (PI)';
hold on

iplot=iplot+1;
plot(xdat_save_pdf_av_low_only{1},y_dat_SW_overall_cont_PD,'b--','linewidth',lwidth);
leg_strs{iplot} = 'Present Day (PD)';


legend(leg_strs,'location','southwest');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
UM_save_plot(gcf,isave_plot_global_cont,savedir,[save_str '_' leg_strs{iplot}] );  


%% Plot the contribution of each cloud fraction to total SW forcing using
%% PD cloud fraction - but only for low-only cloud state and clear state combined
isave_plot_global_cont=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['Forcing_vs_CF'];

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf_av_low_only{iplot});
y_dat_SW_overall_cont_clearlow_PD = ydat_save_pdf_clearlow_only_PD{iplot}.*ydat_save_clearlow_only_PD{iplot}/N;
plot(xdat_save_pdf_clearlow_only_PD{iplot},y_dat_SW_overall_cont_clearlow_PD,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_clearlow_only_PD{iplot},' fraction',' PD');
hold on

legend(leg_strs,'location','southwest');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
UM_save_plot(gcf,isave_plot_global_cont,savedir,[save_str '_' leg_strs{iplot}] );  


% if isave_plot_global_cont==1
%     savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
%     clear opts
%     %        opts.iplot_png=1;
%     opts.iplot_eps=1;
%     saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
% %    close(gcf);
% end 

%% Plot the contribution of each cloud fraction to total SW forcing using
%% PI cloud fractions - but only for low-only cloud state and clear state combined
isave_plot_global_cont=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['Forcing_vs_CF'];

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf_clearlow_only_PI{iplot});
y_dat_SW_overall_cont_clearlow_PI = ydat_save_pdf_clearlow_only_PI{iplot}.*ydat_save_clearlow_only_PI{iplot}/N;
plot(xdat_save_pdf_clearlow_only_PI{iplot},y_dat_SW_overall_cont_clearlow_PI,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_clearlow_only_PI{iplot},' fraction',' PI');
hold on


legend(leg_strs,'location','northwest');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
UM_save_plot(gcf,isave_plot_global_cont,savedir,[save_str '_' leg_strs{iplot}] );  


% if isave_plot_global_cont==1
%     savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
%     clear opts
%     %        opts.iplot_png=1;
%     opts.iplot_eps=1;
%     saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
% %    close(gcf);
% end 

%% PI and PD for the above on the same plot
isave_plot_global_cont=1;
savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['Forcing_vs_CF'];

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0; clear leg_strs

iplot=iplot+1;
plot(xdat_save_pdf_clearlow_only_PI{1},y_dat_SW_overall_cont_clearlow_PI,'b','linewidth',lwidth);
%leg_strs{iplot} = remove_character(xlabelstr_save_pdf_clearlow_only{1},' fraction',' PI');
leg_strs{iplot} = 'Pre-Industrial (PI)';
hold on

iplot=iplot+1;
plot(xdat_save_pdf_clearlow_only_PI{1},y_dat_SW_overall_cont_clearlow_PD,'b--','linewidth',lwidth);
leg_strs{iplot} = 'Present Day (PD)';

% iplot=iplot+1;
% plot(xdat_save_pdf_clearlow_only_PI{1},y_dat_SW_overall_cont_clearlow_PI + y_dat_SW_overall_cont_clearlow_PD,'r','linewidth',lwidth);
% leg_strs{iplot} = 'Sum';

grid on
legend(leg_strs,'location','southwest');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
UM_save_plot(gcf,isave_plot_global_cont,savedir,[save_str '_' leg_strs{iplot}] );  



%% Plot the constribution of each cloud fraction to total SW forcing for PI and PD together
isave_plot_global_cont=0;

% Multiplying the frequency of occurrence by the mean forcing in each bin
% and then normalising
lwidth=4;
fsize=18;

figure
iplot=0;

iplot=iplot+1;
N = sum(ydat_save_pdf_PI{iplot});
plot(xdat_save_pdf_PI{iplot},ydat_save_pdf_PI{iplot}.*ydat_save_PI{iplot}/N,'b','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_PI{iplot},' fraction',' PI');
hold on

iplot=iplot+1;
N = sum(ydat_save_pdf_PI{iplot});
plot(xdat_save_pdf_PI{iplot},ydat_save_pdf_PI{iplot}.*ydat_save_PI{iplot}/N,'r','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_PI{iplot},' fraction',' PI')

iplot=iplot+1;
N = sum(ydat_save_pdf_PI{iplot});
plot(xdat_save_pdf_PI{iplot},ydat_save_pdf_PI{iplot}.*ydat_save_PI{iplot}/N,'k','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf_PI{iplot},' fraction',' PI');

%PD 
iplot=iplot+1;
N = sum(ydat_save_pdf{1});
plot(xdat_save_pdf{1},ydat_save_pdf{1}.*ydat_save{1}/N,'b--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{1},' fraction',' PD');
hold on

iplot=iplot+1;
N = sum(ydat_save_pdf{2});
plot(xdat_save_pdf{2},ydat_save_pdf{2}.*ydat_save{2}/N,'r--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{2},' fraction',' PD')

iplot=iplot+1;
N = sum(ydat_save_pdf{3});
plot(xdat_save_pdf{3},ydat_save_pdf{3}.*ydat_save{3}/N,'k--','linewidth',lwidth);
leg_strs{iplot} = remove_character(xlabelstr_save_pdf{3},' fraction',' PD');

legend(leg_strs,'location','southeast');
xlabel('Cloud fraction');
ylabel('Contribution (W m^{-2})');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-0.05 1.05]);
  


if isave_plot_global_cont==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/Contribution_SW_global'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%    close(gcf);
end 


%% Plots of LWP vs Nd combining oropgraphy and no orography runs

run_set_orog = 'Blending option=3, volcano starting 12 UTC';
run_set_no_orog = 'Blending option=3, volcano starting 12 UTC, no orography';
%run_set_orog = 'Blending option=3, volcano starting 12 UTC, adjusted emission height';
%run_set_no_orog = 'Blending option=3, volcano starting 12 UTC, no orography, adjusted emission height';

plot_type = 'LWP_subgrid PD vs Nd PD';
plot_type = 'CF_subgrid PD vs Nd PD';
plot_type = 'dLWP_subgrid vs dNd'; 
%plot_type = 'dCF_subgrid vs dNd';

SO2_filter_type = 'Volc ON';
CF_filter_str = '';

%Should save the following in the .mat files ideally
Nd_type = 'CASIM';
%xlab_str = ['\DeltaN_{d ' Nd_type '} (cm^{-3})'];
%y_units_str = '(g m^{-2})';
%ylab_str = ['\DeltaLWP ' y_units_str];

    %LON_filter = [-160 -155];
    %LON_filter = [-165 -160];
    %LON_filter = [-170 -165];
    %LON_filter = [-175 -170];
    %LON_filter = [-361 361];    
    
lon_strs={'-160 to -155^{o}E','-165 to -160^{o}E','-170 to -165^{o}E','-175 to -170^{o}E','-361 to 361^{o}E'};
%lon_strs={'-160 to -155^{o}E'};

%-1, 5, 10, 50
LWP_thresh_str = '50 g m^{-2}';
LWP_thresh_str = '10 g m^{-2}';
LWP_thresh_str = '5 g m^{-2}';
LWP_thresh_str = '-1 g m^{-2}';

Nd_type_str='';
%Nd_type_str = '_vs_Rain_Rate_0pt01';
%Nd_type_str = '_vs_Rain_Rate_0pt05';
%Nd_type_str = '_UKCA_Nd';

switch Nd_type_str
    case ''
        exceptions = {'_vs_Rain_Rate','_vs_Rain_Rate_0pt01','_vs_Rain_Rate_0pt05','_UKCA_Nd'};
    otherwise
        exceptions =[];
end

clear leg_str

for ilons=1:length(lon_strs)
    
    lon_str = lon_strs{ilons};
    
    %load in the data
    files=dir([savedir_date run_set_orog ' ' plot_type ' LWP threshold = ' LWP_thresh_str ', SO_2 threshold = 1e-05' SO2_filter_type ', ' CF_filter_str ', ' lon_str ', , minN=*' Nd_type_str '.mat']);
    Hawaii_check_file %script
    loadfile_orog = [savedir_date files(ifile).name];
    dat_orog = load(loadfile_orog);
    
    %files=dir([savedir_date run_set_no_orog ' ' plot_type ' LWP threshold = ' LWP_thresh_str ', SO_2 threshold = 1e-05, ' lon_str ', , minN=*' Nd_type_str '.mat']);
    files=dir([savedir_date run_set_no_orog ' ' plot_type ' LWP threshold = ' LWP_thresh_str ', SO_2 threshold = 1e-05' SO2_filter_type ', ' CF_filter_str ', ' lon_str ', , minN=*' Nd_type_str '.mat']);
    Hawaii_check_file %script
    loadfile_no_orog = [savedir_date files(ifile).name];
    dat_no_orog = load(loadfile_no_orog);
    
    
    figure('color','w');
    %plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
    %mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
    plot(dat_orog.mid_points,dat_orog.yplot,'bo-','linewidth',3); hold on
    leg_str{1} = 'Orography';
    
    xlabel(dat_no_orog.xlab_str);    
    ylabel(dat_no_orog.ylab_str);
    %title(title_str);
    fontsize_figure(gcf,gca,18);
    grid on
    
    plot(dat_no_orog.mid_points,dat_no_orog.yplot,'ro-','linewidth',3);
    leg_str{2} = 'No orography';
    
    legend(leg_str);
    
    title([lon_str ', LWP thresh=' LWP_thresh_str]);
    
    %set(gca,'xlim',[-1000 4000]);
    %set(gca,'ylim',[-600 400]);
            
end
    
%% Lon-height vertical slices / cross sections
var_UM = 'SO2_perkg_lon_height_slice_at_ilat=264_lat=19.50';

um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'SO2_slice_PI_ALL';
UM_load_var_commands

um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'SO2_slice_PD_ALL';
UM_load_var_commands

%Orog off runs

um_case=volc_OFF_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'SO2_slice_PI_ALL2';
UM_load_var_commands

um_case=volc_ON_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'SO2_slice_PD_ALL2';
UM_load_var_commands

%

var_UM = 'LWC_lon_height_slice_at_ilat=264_lat=19.50';

um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'LWC_slice_PI_ALL';
UM_load_var_commands

um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'LWC_slice_PD_ALL';
UM_load_var_commands

um_case=volc_OFF_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'LWC_slice_PI_ALL2';
UM_load_var_commands

um_case=volc_ON_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'LWC_slice_PD_ALL2';
UM_load_var_commands


%

var_UM = 'potential_temperature_lon_height_slice_at_ilat=264_lat=19.50';

um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'potemp_slice_PI_ALL';
UM_load_var_commands

um_case=um_case_PD; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'potemp_slice_PD_ALL';
UM_load_var_commands

%Orog off
um_case=volc_OFF_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'potemp_slice_PI_ALL2';
UM_load_var_commands

um_case=volc_ON_no_orog; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'potemp_slice_PD_ALL2';
UM_load_var_commands


%Orography height
var_UM = 'Orography_lon_height_slice_at_ilat=264_lat=19.50';

um_case='u-ch764'; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
opts.lat_var='height'; opts.lon_var='Longitude';
var_name_out = 'orog_PD_ALL'; icheck_time=0;
UM_load_var_commands

   
%% Plot slices - sort the height, orography, etc.


z = dat_global.gcm_Plat2D_edges_UM(:,1); %lat is used for height for the slices
z(1)=0;
%zedges = [0; 0.5* (z(1:end-1) + z(2:end) )];
x = dat_global.gcm_Plon2D_edges_UM(1,:);

max_z = 8e3;
iz = find(z<=max_z); iz=iz(end);
zplot = z(1:iz)/1000; zplot = zplot';
[X,Z] = meshgrid(x,zplot);
%Can do this to add the orography:-
lon_mids = squeeze(dat_global.gcm_Plon2D_UM(1,:));
orog_mids = squeeze(orog_PD_ALL(1,:))'/1e3;
sample_pts = squeeze(dat_global.gcm_Plon2D_edges_UM(1,:));
orog = interp1(lon_mids,orog_mids,sample_pts,'linear','extrap');
orog2d = repmat(orog,[size(Z,1) 1]);
Z = Z + orog2d;
%Check whether this is the correct way to do the heights by calculating
%properly.


%% SO2 slice

it=1;
%it=7;
%it=96;
% SO2
dat = squeeze(SO2_slice_PD_ALL(1:iz-1,:,it));
figure('color','w');
%dpcolor(x,zplot,dat);
dpcolor(X,Z,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
%set(gca,'xlim',[-157 -153]);
set(gca,'ylim',[0 zplot(end)]);
title(['Volcano ON, orography ON it=' num2str(it)]);
caxis([0 2e-6]);

dat = squeeze(SO2_slice_PD_ALL2(1:iz-1,:,it));
figure('color','w');
dpcolor(x,zplot,dat);
%dpcolor(X,Z,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
%set(gca,'xlim',[-157 -153]);
set(gca,'ylim',[0 zplot(end)]);
title(['Volcano ON, orography OFF it=' num2str(it)]);
caxis([0 2e-6]);


%% LWC

it=1;
it=7;
%it=96;

dat = squeeze(LWC_slice_PD_ALL(1:iz-1,:,it));
figure('color','w');
dpcolor(X,Z,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
set(gca,'xlim',[-161 -153]);
set(gca,'ylim',[0 zplot(end)]);
title(['LWC; Volcano ON, orography ON it=' num2str(it)]);
caxis([0 1e-3]);

dat = squeeze(LWC_slice_PD_ALL2(1:iz-1,:,it));
figure('color','w');
dpcolor(x,zplot,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
set(gca,'xlim',[-161 -153]);
set(gca,'ylim',[0 zplot(end)]);
title(['LWC; Volcano ON, orography OFF it=' num2str(it)]);
caxis([0 1e-3]);

%% LWC time mean

dat = meanNoNan(LWC_slice_PD_ALL(1:iz-1,:,:),3);
figure('color','w');
dpcolor(X,Z,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
set(gca,'xlim',[-161 -153]);
set(gca,'ylim',[0 zplot(end)]);
title('LWC; Volcano ON, orography ON, time mean');
caxis([0 1.5e-4]);

dat = meanNoNan(LWC_slice_PD_ALL2(1:iz-1,:,:),3);
figure('color','w');
dpcolor(x,zplot,dat);
%dpcolor(X,Z,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
set(gca,'xlim',[-161 -153]);
set(gca,'ylim',[0 zplot(end)]);
title('LWC; Volcano ON, orography OFF,  time mean');
caxis([0 1.5e-4]);


%% Potential temperature slice


it=1;
it=7;
%it=96;
% SO2
dat = squeeze(potemp_slice_PD_ALL(1:iz-1,:,it));
figure('color','w');
%dpcolor(x,zplot,dat);
dpcolor(X,Z,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
set(gca,'xlim',[-161 -153]);
set(gca,'ylim',[0 zplot(end)]);
title(['Volcano ON, orography ON it=' num2str(it)]);
caxis([295 330]);

dat = squeeze(potemp_slice_PD_ALL2(1:iz-1,:,it));
figure('color','w');
dpcolor(x,zplot,dat);
%dpcolor(X,Z,dat); 
shading flat; colorbar
xlabel('Longitude (degrees)');
ylabel('Height (km)');
set(gca,'xlim',[-161 -153]);
set(gca,'ylim',[0 zplot(end)]);
title(['Volcano ON, orography OFF it=' num2str(it)]);
caxis([295 330]);



%% Single line plots across the plume

lat_plume_cs = [28 8]; dlat_cs = diff(lat_plume_cs);
lon_plume_cs = [-165 -160]; dlon_cs = diff(lon_plume_cs);

npts = 300;
lat_plume_cs_single = [lat_plume_cs(1):dlat_cs/npts:lat_plume_cs(2)];
lon_plume_cs_single = [lon_plume_cs(1):dlon_cs/npts:lon_plume_cs(2)];

cs_Nd_PD = griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,Nd_PD_map,lon_plume_cs_single(:),lat_plume_cs_single(:));
cs_SO2_PD = griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,SO2_PD,lon_plume_cs_single(:),lat_plume_cs_single(:));
cs_LWP_PD = griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,mean_LWP_PD,lon_plume_cs_single(:),lat_plume_cs_single(:));
cs_LWP_PI = griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,mean_LWP_PI,lon_plume_cs_single(:),lat_plume_cs_single(:));

figure
plot(lat_plume_cs_single,cs_Nd_PD);

figure
plot(lat_plume_cs_single,cs_SO2_PD);

figure
hold on
plot(lat_plume_cs_single,cs_LWP_PI,'r');
plot(lat_plume_cs_single,cs_LWP_PD);

%% Multiple interpolated lines across the plume. It's actually probably just as useful to look at the maps
% of the time-mean LWP and to try and pick out a plume by eye. But can
% do this using only data when SO2>threshold to allow for the time
% evolution of the plume (accounting for the fact that the plume has not
% reached the points further downwind early on in the sim).
% For satellite data of an ongoing plume perhaps it is less of a concern.

%This defines points perpdenicular to the wind/plume direction - i.e.,
%across the plume.
lat_plume_cs = [28 8]; dlat_cs = diff(lat_plume_cs);
lon_plume_cs = [-165 -160]; dlon_cs = diff(lon_plume_cs);

%Multiple points along the across-plume line - these will be the start
%points of high resolution along-plume lines for averaging
npts = 50;
lat_plume_cs2 = [lat_plume_cs(1):dlat_cs/npts:lat_plume_cs(2)];
lon_plume_cs2 = [lon_plume_cs(1):dlon_cs/npts:lon_plume_cs(2)];
%x and y distances moved by the line in physical distance rather than lat
%lon - not sure that this is necessary?
x = distlatlon(lat_plume_cs(1),lon_plume_cs(1),lat_plume_cs(1),lon_plume_cs(2));
y = distlatlon(lat_plume_cs(2),lon_plume_cs(2),lat_plume_cs(1),lon_plume_cs(2));
%Angle relative to north (anti-clockwise) :-
th_across = atan(x./y);
%Angle in lat/lon space (alternative, simpler method).
th_across2 = atan( (lon_plume_cs(2) - lon_plume_cs(1)) ./ (lat_plume_cs(2) - lat_plume_cs(1)) );


%Make a matrix of all the along-plume lines combined
D = 2000; %Distance to move downwind
x2 = D.*cos(th_across); %dx change for the downwind vector
y2 = D.*sin(th_across);
dlat = y2 / 111.3195; %using 1 deg as 6378.140*2*pi/360 - same earth radius as in distlatlon

%Simpler method. Convert the distance downwind (D) into an approx distance in lat/lon space :-
dlon2 = D/111.3195.*cos(th_across2); %dx change for the downwind vector
dlat2 = -D/111.3195.*sin(th_across2);


npts=100;
lats_ALL = NaN*ones([length(lat_plume_cs2) npts]);
lons_ALL = NaN*ones([length(lat_plume_cs2) npts]);
for i=1:length(lat_plume_cs2)
    %cs_LWP_PD_multi(ilon,:) = griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,mean_LWP_PD,lon_plume_cs2(:)+dlon_multi*ilon,lat_plume_cs2(:));
    
    latA = lat_plume_cs2(i);
    latB = lat_plume_cs2(i) - dlat;
    dlon = distlatlon_inverse(x2,latA,latB);    
        
    %dlon=dlon2; dlat=dlat2; %if testing simpler method.
    lons = [lon_plume_cs2(i):-dlon/npts:lon_plume_cs2(i)-dlon];
    lats = [lat_plume_cs2(i):-dlat/npts:lat_plume_cs2(i)-dlat];
    lats_ALL(i,1:npts) = lats(1:npts);
    lons_ALL(i,1:npts) = lons(1:npts);        
end

cs_LWP_PD_multi = griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,mean_LWP_PD,lons_ALL,lats_ALL);
cs_LWP_PD_multi_mean = meanNoNan(cs_LWP_PD_multi,2);

figure
plot(lat_plume_cs2,cs_LWP_PD_multi_mean);

cs_Nd_PD_multi = griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,Nd_PD_map,lons_ALL,lats_ALL);
[cs_Nd_PD_multi_mean,N,s] = meanNoNan(cs_Nd_PD_multi,2);
cs_Nd_PD_multi_mean(N<60)=NaN;

figure
plot(lat_plume_cs2,cs_Nd_PD_multi_mean);


%single plume
figure
plot(lat_plume_cs_single,cs_Nd_PD);



%% time mean dLWP only using times when the SO2 is above the threshold

%     mean_LWP_PI = meanNoNan(LWP_PI_ALL,3);
%     mean_LWP_PD = meanNoNan(LWP_PD_ALL,3); 
%     mean_LWP_PI2 = meanNoNan(LWP_PI_ALL2,3);
%     mean_LWP_PD2 = meanNoNan(LWP_PD_ALL2,3); 
    
%Need to also do this for no-orog run
    thresh_SO2 = 1e-5; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);
    thresh_SO2 = 1e-4; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);

    lwp_SO2_PI = LWP_PI_ALL; lwp_SO2_PI(inan) = NaN;
    lwp_SO2_PD = LWP_PD_ALL; lwp_SO2_PD(inan) = NaN;
    lwp_SO2_PI2 = LWP_PI_ALL2; lwp_SO2_PI2(inan) = NaN;
    lwp_SO2_PD2 = LWP_PD_ALL2; lwp_SO2_PD2(inan) = NaN;
    
    [mean_lwp_SO2_PI,N_mean_lwp_SO2_PI std_dev] = meanNoNan(lwp_SO2_PI,3);
    [mean_lwp_SO2_PD,N_mean_lwp_SO2_PD std_dev] = meanNoNan(lwp_SO2_PD,3);
    [mean_lwp_SO2_PI2,N_mean_lwp_SO2_PI2 std_dev] = meanNoNan(lwp_SO2_PI2,3);
    [mean_lwp_SO2_PD2,N_mean_lwp_SO2_PD2, std_dev] = meanNoNan(lwp_SO2_PD2,3);
    
    thresh_N = 0; %Overwritten below
    icoarse_grain=1;
    M_coarse_grain=10; N_coarse_grain=10;    
    var_UM = '';
    dat_modis = mean_lwp_SO2_PD - mean_lwp_SO2_PI;
    thresh_N = 10; dat_modis(N_mean_lwp_SO2_PI<thresh_N) = NaN;
    subtitle_str = ['Time mean \DeltaLWP (g m^{-2}) for SO_2>' num2str(thresh_SO2,'%1.0e') ', N>=' num2str(thresh_N)];
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]);
    icoarse_grain=0;        
    
    thresh_N = 0; %Overwritten below
    icoarse_grain=1;
    M_coarse_grain=10; N_coarse_grain=10;   
    var_UM = '';
    dat_modis = mean_lwp_SO2_PD2 - mean_lwp_SO2_PI2;
    thresh_N = 10; dat_modis(N_mean_lwp_SO2_PI<thresh_N) = NaN;
    subtitle_str = ['Time mean \DeltaLWP (g m^{-2}) for SO_2>' num2str(thresh_SO2,'%1.0e') ' NO OROG, N>=' num2str(thresh_N)];
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]);
    icoarse_grain=0;
    
    icoarse_grain=1;
    M_coarse_grain=10; N_coarse_grain=10;
    subtitle_str = ['Npoints for SO_2>' num2str(thresh_SO2,'%1.0e')];
    var_UM = '';
    dat_modis = N_mean_lwp_SO2_PI;
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 100]);
    icoarse_grain=0; 
    
%%  Time mean dLWP, also filtering for LWP. Used for the PDFs in the next step

%if exist('f0_orig')
    Hawaii_calc_LWPic_using_subgrid_CF %calculates the LWPic using the subgrid CF
%end

%also filter out LWP<50 and dNd<50.
%PI_or_PD_thresh_LWP = 'both';
%PI_or_PD_thresh_LWP = 'PI';

isingle_time=0; it=33; 


%vars_multi = {'LWPic_using_CF_subgrid_vert_int','CF_subgrid_vert_int','Nd_using_CF_subgrid_vert_int','LWP','CF tau thresh'...
%    'CF','LWPic'};
vars_multi = {'LWPic_using_CF_subgrid_vert_int specific time'};
vars_multi = {'LWPic_using_CF_subgrid_vert_int'};
%vars_multi = {'CF_subgrid_vert_int'};
%vars_multi = {'Nd_using_CF_subgrid_vert_int'};
%vars_multi = {'LWP'};
%vars_multi = {'Cloud base rain rates'};
%vars_multi = {'RWP'};
%vars_multi = {'Liquid Cloud Top Height'};
%vars_multi = {'Total Ice Water Path'};
%vars_multi = {'Ice Water Path Section 30'};
%vars_multi = {'total_max_random_cloud_amount'};

clear  mean_overall N_overall std_dev_overall mean_overall_prc

for ivar=1:length(vars_multi)
    
    var_mean = vars_multi{ivar};

%var_mean = 'LWP';
%var_mean = 'LWPic';
%var_mean = 'Nd';
%var_mean = 'CF';
%var_mean = 'CF tau thresh';
%var_mean = 'LWP tau thresh';
%var_mean = 'CF_subgrid_vert_int';
%var_mean = 'LWPic_using_CF_subgrid_vert_int';
%var_mean = 'Nd_using_CF_subgrid_vert_int';
%var_mean = 'SW_estimated';
%var_mean = 'SW_actual';
%var_mean = 'ERF_CF';
%var_mean = 'ERF_LWPic';
%var_mean = 'ERF_Nd';
%var_mean = 'ERF_tot Estimated';
%var_mean = 'ERF_sum Estimated';
%var_mean = 'ERF_aci Actual';
%var_mean = 'ERF_ari Actual';


plot_type = 'delta';
%plot_type = 'PD';
%plot_type = 'PI';

PI_or_PD_thresh_LWP = 'none'; %Prob best setting none as the default and then using this only for in-cloud variables (LWP, Nd)

dSO2 = SO2_col_PD_ALL - SO2_col_PI_ALL;

thresh_SO2(1) = 1e-5; thresh_SO2(2) = 1e15; inan = find(dSO2 < thresh_SO2(1) | dSO2 >= thresh_SO2(2));
%thresh_SO2(1) = -1e-5; thresh_SO2(2) = 1e-5; inan = find(dSO2 < thresh_SO2(1) | dSO2 >= thresh_SO2(2));
%thresh_SO2(1) = -1e-5; thresh_SO2(2) = 1e-5; inan = find(SO2_col_PD_ALL < thresh_SO2(1) | SO2_col_PD_ALL >= thresh_SO2(2));
%thresh_SO2(1) = 1e-6; thresh_SO2(2) = 1e-5; inan = find(SO2_col_PD_ALL < thresh_SO2(1) | SO2_col_PD_ALL >= thresh_SO2(2));
%thresh_SO2(1) = 1e-5; thresh_SO2(2) = 1e15; inan = find(dSO2 < thresh_SO2(1) | dSO2 >= thresh_SO2(2));
%thresh_SO2(1) = -1e-5; thresh_SO2(2) = 1e-5; inan = find(dSO2 < thresh_SO2(1) | dSO2 >= thresh_SO2(2));
%thresh_SO2 = 1e-5; inan = find(SO2_col_PD_ALL < thresh_SO2);
%thresh_SO2 = 1e-5; inan = find(dSO2 < thresh_SO2);

%thresh_SO2 = 1e-4; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);
%thresh_SO2 = -1e9; inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);



switch plot_type
    case 'delta'
        i_delta=1;
    otherwise
        i_delta=0;
end



CTT=278; P=850e2; k=0.8;
[tau_PI_ALL,H_PI]=MODIS_tau_func_N_LWP(LWP_PI_ALL*1e-3,Nd_PI_ALL,CTT,k,P); %LWP needs to be in kg/m2 and Nd in per m3
[tau_PD_ALL,H_PD]=MODIS_tau_func_N_LWP(LWP_PD_ALL*1e-3,Nd_PD_ALL,CTT,k,P);


clims = [-150 150];
 
        

switch var_mean
    case 'total_max_random_cloud_amount'
        dat_SO2_PI =  total_max_random_cloud_amount_PI_ALL; dat_SO2_PI(inan) = NaN; %g/m2
        dat_SO2_PD =  total_max_random_cloud_amount_PD_ALL; dat_SO2_PD(inan) = NaN;
                
        %dat_SO2_PI = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt05; dat_SO2_PI(inan) = NaN; %mm hr^-1
        %dat_SO2_PD = 3600*RainRate_0pt05_PD_ALL; dat_SO2_PD(inan) = NaN;
        
        var_str = 'Total Cloud Fraction (max random overlap)';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];
        
        %min_val = 1e-8;
        %dat_SO2_PI(dat_SO2_PI<min_val)=min_val;
        %dat_SO2_PD(dat_SO2_PD<min_val)=min_val;
        
    case 'Ice Water Path Section 30'
        dat_SO2_PI = IWP_sec30_PI_ALL *1e3; dat_SO2_PI(inan) = NaN; %g/m2
        dat_SO2_PD = IWP_sec30_PD_ALL*1e3; dat_SO2_PD(inan) = NaN;
                
        %dat_SO2_PI = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt05; dat_SO2_PI(inan) = NaN; %mm hr^-1
        %dat_SO2_PD = 3600*RainRate_0pt05_PD_ALL; dat_SO2_PD(inan) = NaN;
        
        var_str = 'Ice water path section 30 (g m^{-2})';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];
        
        min_val = 1e-8;
        dat_SO2_PI(dat_SO2_PI<min_val)=min_val;
        dat_SO2_PD(dat_SO2_PD<min_val)=min_val;
        
    case 'Total Ice Water Path'
        dat_SO2_PI = TIWP_total_column_to_zdomain_top_PI_ALL*1e3; dat_SO2_PI(inan) = NaN; %g/m2
        dat_SO2_PD = TIWP_total_column_to_zdomain_top_PD_ALL*1e3; dat_SO2_PD(inan) = NaN;
                
        %dat_SO2_PI = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt05; dat_SO2_PI(inan) = NaN; %mm hr^-1
        %dat_SO2_PD = 3600*RainRate_0pt05_PD_ALL; dat_SO2_PD(inan) = NaN;
        
        var_str = 'Ice water path (g m^{-2})';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80]; 
        
        min_val = 1e-8;
        dat_SO2_PI(dat_SO2_PI<min_val)=min_val;
        dat_SO2_PD(dat_SO2_PD<min_val)=min_val;
    
    case 'Liquid Cloud Top Height'
        dat_SO2_PI = Max_Cloud_Height_in_cloud_LWC_PI_ALL/1e3; dat_SO2_PI(inan) = NaN; %km
        dat_SO2_PD = Max_Cloud_Height_in_cloud_LWC_PD_ALL/1e3; dat_SO2_PD(inan) = NaN;
                
        %dat_SO2_PI = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt05; dat_SO2_PI(inan) = NaN; %mm hr^-1
        %dat_SO2_PD = 3600*RainRate_0pt05_PD_ALL; dat_SO2_PD(inan) = NaN;
        
        var_str = 'Liquid Cloud Top Height (km)';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];  
        
    case 'RWP'
        dat_SO2_PI = 1e3*RWP_total_column_to_zdomain_top_PI_ALL; dat_SO2_PI(inan) = NaN; %g m^-2
        dat_SO2_PD = 1e3*RWP_total_column_to_zdomain_top_PD_ALL; dat_SO2_PD(inan) = NaN;
                
        %dat_SO2_PI = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt05; dat_SO2_PI(inan) = NaN; %mm hr^-1
        %dat_SO2_PD = 3600*RainRate_0pt05_PD_ALL; dat_SO2_PD(inan) = NaN;
        
        var_str = 'Rain Water Path (g m^{-2})';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];  
        
     case 'Cloud base rain rates'
        %dat_SO2_PI = 3600*RainRate_0pt01_PI_ALL; dat_SO2_PI(inan) = NaN; %mm hr^-1
        %dat_SO2_PD = 3600*RainRate_0pt01_PD_ALL; dat_SO2_PD(inan) = NaN;
                
        %dat_SO2_PI = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt01_PI_ALL; dat_SO2_PI(inan) = NaN; %mm hr^-1
        %dat_SO2_PD = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt01_PD_ALL; dat_SO2_PD(inan) = NaN;
        
        dat_SO2_PI = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt01_PI_ALL; dat_SO2_PI(inan) = NaN; %mm hr^-1
        dat_SO2_PD = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt01_PD_ALL; dat_SO2_PD(inan) = NaN;
        

        
        var_str = 'Cloud base rain rate (>0.05 g kg^{-1} LWC threshold)';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];    
        
        min_val = 1e-6;
        dat_SO2_PI(dat_SO2_PI<min_val)=min_val;
        dat_SO2_PD(dat_SO2_PD<min_val)=min_val;
        
     case 'ERF_ari Actual'
        dat_SO2_PI = ERF_ARI_TOA_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = ERF_ARI_TOA_ALL; dat_SO2_PD(inan) = NaN;
                
        var_str = 'ERF_{ari} Actual';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];  
        
     case 'ERF_aci Actual'
        dat_SO2_PI = ERF_ACI_TOA_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = ERF_ACI_TOA_ALL; dat_SO2_PD(inan) = NaN;
        
        
        var_str = 'ERF_{aci} Actual';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];         
         
     case 'ERF_sum Estimated'
        %dat_SO2_PI = SW_estimated_forcing_sum_TOA_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = SW_estimated_forcing_sum_PD_TOA_ALL; dat_SO2_PD(inan) = NaN;
        dat_SO2_PI = SW_estimated_forcing_sum_PIPD_TOA_ALL; dat_SO2_PI(inan) = NaN;
        
        var_str = 'ERF_{sum} Estimated';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];
        
    case 'ERF_tot Estimated'
        dat_SO2_PI = SW_estimated_forcing_tot_TOA; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = SW_estimated_forcing_tot_TOA; dat_SO2_PD(inan) = NaN;
        var_str = 'ERF_{tot} Estimated';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];
        
    case 'ERF_CF'
        dat_SO2_PI = SW_estimated_forcing_cf_PIPD_TOA_ALL; dat_SO2_PI(inan) = NaN;
        at_SO2_PD = SW_estimated_forcing_cf_PIPD_TOA_ALL; dat_SO2_PD(inan) = NaN;  
        
        %dat_SO2_PI = SW_estimated_forcing_cf_TOA; dat_SO2_PI(inan) = NaN;
        %dat_SO2_PD = SW_estimated_forcing_cf_PD_TOA; dat_SO2_PD(inan) = NaN;
        
        var_str = 'ERF_{fc}';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];
        
    case 'ERF_LWPic';        
        
        dat_SO2_PI = SW_estimated_forcing_lwp_PIPD_TOA_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = SW_estimated_forcing_lwp_PIPD_TOA_ALL; dat_SO2_PD(inan) = NaN;  
        
        %dat_SO2_PI = SW_estimated_forcing_lwp_TOA; dat_SO2_PI(inan) = NaN;
        %dat_SO2_PD = SW_estimated_forcing_lwp_PD_TOA; dat_SO2_PD(inan) = NaN;
        
        var_str = 'ERF_{LWPic}';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];
        
    case 'ERF_Nd';
        dat_SO2_PI = SW_estimated_forcing_Nd_PIPD_TOA_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = SW_estimated_forcing_Nd_PIPD_TOA_ALL; dat_SO2_PD(inan) = NaN;
        %dat_SO2_PI = SW_estimated_forcing_Nd_TOA; dat_SO2_PI(inan) = NaN;
        %dat_SO2_PD = SW_estimated_forcing_Nd_PD_TOA; dat_SO2_PD(inan) = NaN;
                
        var_str = 'ERF_{Nd}';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];
        
     case 'SW_estimated'
        dat_SO2_PI = SW_PI_TOA ; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = SW_PD_TOA; dat_SO2_PD(inan) = NaN;
        var_str = 'SW_{TOA} Estimated';
        if i_delta==1
            var_str = ['\Delta' var_str];        
        end
        PI_or_PD_thresh_LWP = 'none';  
        clims = [-80 80];
        
    case 'SW_actual'
        dat_SO2_PI = SW_up_TOA_PI_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = SW_up_TOA_PD_ALL; dat_SO2_PD(inan) = NaN;
        var_str = 'SW_{TOA} Actual';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';
        clims = [-80 80];    
        
        
    case 'LWP'
        dat_SO2_PI = LWP_PI_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = LWP_PD_ALL; dat_SO2_PD(inan) = NaN;
        if i_delta==1
            var_str = '\DeltaLWP (g m^{-2})';
        else
            var_str = 'LWP (g m^{-2})';
        end
        %PI_or_PD_thresh_LWP = 'both';
        PI_or_PD_thresh_LWP = 'none';        
        
    case 'LWPic'
        dat_SO2_PI = LWP_PI_ALL; dat_SO2_PI(inan) = NaN; %Is this wrong?
        dat_SO2_PD = LWP_PD_ALL; dat_SO2_PD(inan) = NaN;
        if i_delta==1
            var_str = '\DeltaLWP (g m^{-2})';
        else
            var_str = 'LWP (g m^{-2})';
        end
        PI_or_PD_thresh_LWP = 'both';
        %PI_or_PD_thresh_LWP = 'none';
        
    case 'LWP tau thresh'
        dat_SO2_PI = LWP_PI_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = LWP_PD_ALL; dat_SO2_PD(inan) = NaN;
        if i_delta==1
            var_str = '\DeltaLWP (g m^{-2})';
        else
            var_str = 'LWP (g m^{-2})';
        end
        
        PI_or_PD_thresh_LWP = 'both using tau';
        
    case 'Nd'
        dat_SO2_PI = Nd_PI_ALL/1e6; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = Nd_PD_ALL/1e6; dat_SO2_PD(inan) = NaN;
        if i_delta==1
            var_str = '\DeltaN_d (cm^{-3})';
        else
            var_str = 'N_d (cm^{-3})';
        end
        PI_or_PD_thresh_LWP = 'both';
        
    case 'CF'
        thresh_lwp_CF = 5;
        thresh_lwp_CF = 1; %LWP=1 corresponds to tau=0.32 for Nd=100 per cc.
        
        dat_SO2_PI = zeros(size(LWP_PI_ALL)); %dat_SO2_PI(inan) = NaN;
        dat_SO2_PI(isnan(LWP_PI_ALL)) = NaN;
        dat_SO2_PI(LWP_PI_ALL>thresh_lwp_CF) = 1;
        dat_SO2_PI(inan) = NaN;
        
        dat_SO2_PD = zeros(size(LWP_PD_ALL)); %dat_SO2_PD(inan) = NaN;
        dat_SO2_PD(isnan(LWP_PD_ALL)) = NaN;
        dat_SO2_PD(LWP_PD_ALL>thresh_lwp_CF) = 1;
        dat_SO2_PD(inan) = NaN;
        
        if i_delta==1
            var_str = '\DeltaCF';    
        else
            var_str = 'CF';
        end
        
        PI_or_PD_thresh_LWP = 'none';
        
    case 'CF tau thresh'
        %Estimate COD using LWP and Nd.
        
        thresh_lwp_tau = 0.3;
        
              
        
        dat_SO2_PI = zeros(size(LWP_PI_ALL)); %dat_SO2_PI(inan) = NaN;
        dat_SO2_PI(isnan(LWP_PI_ALL)) = NaN;
        dat_SO2_PI(tau_PI_ALL>thresh_lwp_tau) = 1;
        dat_SO2_PI(inan) = NaN;
        
        dat_SO2_PD = zeros(size(LWP_PD_ALL)); %dat_SO2_PD(inan) = NaN;
        dat_SO2_PD(isnan(LWP_PD_ALL)) = NaN;
        dat_SO2_PD(tau_PD_ALL>thresh_lwp_tau) = 1;
        dat_SO2_PD(inan) = NaN;
        
        cf_str = ['CF_{tau>' num2str(thresh_lwp_tau) '}'];
        if i_delta==1
            var_str = ['\Delta' cf_str];
        else
            var_str = cf_str;
        end   
                    
    case 'CF_subgrid_vert_int' 
        dat_SO2_PI = CF_lwc_weighted_total_column_to_zdomain_top_PI_ALL; dat_SO2_PI(inan) = NaN;
        dat_SO2_PD = CF_lwc_weighted_total_column_to_zdomain_top_PD_ALL; dat_SO2_PD(inan) = NaN;
        if i_delta==1
            var_str = '\DeltaCF_{vi}';
        else
            var_str = 'CF_{vi}';
        end
        PI_or_PD_thresh_LWP = 'none';   
        
    case 'LWPic_using_CF_subgrid_vert_int'
        dat_SO2_PI = W0_orig; dat_SO2_PI(inan) = NaN; %W0_orig calculated in Hawaii_calc_LWPic_using_subgrid_CF
        dat_SO2_PD = W1_orig; dat_SO2_PD(inan) = NaN;
        var_str='LWP_{ic sub-grid}';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';      
        
case 'LWPic_using_CF_subgrid_vert_int specific time'
        dat_SO2_PI = W0_orig; dat_SO2_PI(inan) = NaN; %W0_orig calculated in Hawaii_calc_LWPic_using_subgrid_CF
        dat_SO2_PD = W1_orig; dat_SO2_PD(inan) = NaN;
        isingle_time=1; it=33;
        var_str=['LWP_{ic sub-grid}, ' datestr(time_out(it))];
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none';         
        
    case 'Nd_using_CF_subgrid_vert_int'
        dat_SO2_PI = N0_orig; dat_SO2_PI(inan) = NaN; %W0_orig calculated in Hawaii_calc_LWPic_using_subgrid_CF
        dat_SO2_PD = N1_orig; dat_SO2_PD(inan) = NaN;
        var_str='Nd_{CF>0.01}';
        if i_delta==1
            var_str = ['\Delta' var_str];
        end
        PI_or_PD_thresh_LWP = 'none'; 
        
   
end
    


switch PI_or_PD_thresh_LWP
    case 'both'
        thresh_LWP = 50; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
        %thresh_LWP = 10; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
        thresh_LWP = 5; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
        thresh_LWP = 1; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
        %thresh_LWP = -1; inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
        %var_test(inan)=NaN; var_filter(inan)=NaN; %Remove the filter values too to ensure even bin
        %sampling when using percentiles
        
        thresh_str = ['LWP>' num2str(thresh_LWP)];
        
    case 'PI'
        %Don't think this is a fair comparison since only the high PI
        %values are left, but all the low PD ones remain.
        thresh_LWP = 10; inan = find(LWP_PI_ALL<thresh_LWP);
        
        thresh_str = ['LWP>' num2str(thresh_LWP)];
        
    case 'both using tau'        
        thresh_LWP = 0.3; inan = find(tau_PI_ALL<thresh_LWP | tau_PD_ALL<thresh_LWP);   
        
        thresh_str = ['\tau>' num2str(thresh_LWP)];
        
    otherwise
        inan=[];
        
end
        
 dat_SO2_Nd_PI = dat_SO2_PI; dat_SO2_Nd_PI(inan) = NaN;
 dat_SO2_Nd_PD = dat_SO2_PD; dat_SO2_Nd_PD(inan) = NaN;
 %dat_SO2_PI2(inan) = NaN;
 %dat_SO2_PD2(inan) = NaN;

 % Filter by dNd too
 dNd = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
 thresh_dNd = -1e9; inan = find(dNd < thresh_dNd);
 %thresh_dNd = 10; inan = find(dNd < thresh_dNd);
 %thresh_dNd = 50; inan = find(dNd < thresh_dNd);
 
 dat_SO2_Nd_PI(inan) = NaN;
 dat_SO2_Nd_PD(inan) = NaN;
 %dat_SO2_PI2(inan) = NaN;
 %dat_SO2_PD2(inan) = NaN;
 
    
    %[mean_dat_SO2_dNd_PI,N_mean_dat_SO2_dNd_PI std_dev] = meanNoNan(dat_SO2_Nd_PI,3);
    %[mean_dat_SO2_dNd_PD,N_mean_dat_SO2_dNd_PD std_dev] = meanNoNan(dat_SO2_Nd_PD,3);
    %[mean_dat_SO2_dNd_PI2,N_mean_dat_SO2_dNd_PI2 std_dev] = meanNoNan(dat_SO2_PI2,3);
    %[mean_dat_SO2_dNd_PD2,N_mean_dat_SO2_dNd_PD2, std_dev] = meanNoNan(dat_SO2_PD2,3); 
    
    switch plot_type
        case 'delta'            
            d_dat = dat_SO2_Nd_PD - dat_SO2_Nd_PI;
        case 'PI'
            d_dat = dat_SO2_Nd_PI;
        case 'PD'
            d_dat = dat_SO2_Nd_PD;
    end
    
    if isingle_time==1
        d_dat = d_dat(:,:,it);
    end
    
    [mean_d_dat,N_mean_d_dat std_dev] = meanNoNan(d_dat,3);
    mean_d_dat_prc = 100 * mean_d_dat ./ meanNoNan(dat_SO2_Nd_PI,3);
    
    [mean_overall(ivar),N_overall(ivar) std_dev_overall(ivar)] = meanNoNan(d_dat(:),1)
    
    mean_overall_prc(ivar) = 100 * mean_overall(ivar) / meanNoNan(dat_SO2_Nd_PI(:),1)
    
    
    
    
 thresh_N = 0; %Overwritten below
 icoarse_grain=1;
 M_coarse_grain=10; N_coarse_grain=10;
 var_UM = '';
 %dat_modis = mean_dat_SO2_dNd_PD - mean_dat_SO2_dNd_PI;
 %thresh_N = 5; dat_modis(N_mean_dat_SO2_dNd_PI<thresh_N) = NaN;
 dat_modis = mean_d_dat;
 %thresh_N = 5; dat_modis(N_mean_d_dat<thresh_N) = NaN;
 subtitle_str = ['Time mean ' var_str ' for SO_2>' num2str(thresh_SO2,'%1.0e') ', ' thresh_str ...
     ', \DeltaN_d>' num2str(thresh_dNd) ', N>=' num2str(thresh_N)];
 figure
 UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
 lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
 caxis(clims);
 icoarse_grain=0;
 
 thresh_N = 0; %Overwritten below
 icoarse_grain=1;
 M_coarse_grain=10; N_coarse_grain=10;
 var_UM = '';
 %dat_modis = mean_dat_SO2_dNd_PD - mean_dat_SO2_dNd_PI;
 %thresh_N = 5; dat_modis(N_mean_dat_SO2_dNd_PI<thresh_N) = NaN;
 dat_modis = mean_d_dat_prc;
 %thresh_N = 5; dat_modis(N_mean_d_dat<thresh_N) = NaN;
 subtitle_str = ['% Time mean ' var_str ' for SO_2>' num2str(thresh_SO2,'%1.0e') ', ' thresh_str ...
     ', \DeltaN_d>' num2str(thresh_dNd) ', N>=' num2str(thresh_N)];
 figure
 UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
 lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
 caxis(clims);
 icoarse_grain=0;
 
end
 
 %% Make a 1D PDF of the delta data from the above
 PDF_type = 'delta time mean';
 PDF_type = 'PI and PD';
 
 if ~exist('thresh_LWP')
     thresh_LWP = -1;
 end
 
 logbin_norm_driver=0; %set to one if using log bins and setting the xscale to log so that the area
 %Under the PDFs in log space will be the no. datapoints.
 
 switch PDF_type
     case {'delta all data','delta time mean'}
         switch var_mean
             case {'LWP','LWPic_using_CF_subgrid_vert_int'}
                 Ybins_DRIVER = [-500:10:500];
                 ylabelstr='\DeltaLWP_{cloudy sky} (g m^{-2})';
             case 'Nd'
                 Ybins_DRIVER = [-500:10:500];
                 ylabelstr='\DeltaNd (cm^{-3})';
             case 'CF'
                 Ybins_DRIVER = [-1:0.01:1];
                 ylabelstr='\DeltaCF';              
         end
         
         switch PDF_type
             case 'delta all data'                 
                 Y_driver = d_dat(:);
                 
             case 'delta time mean'
                 Y_driver = mean_d_dat(:);   
                 %Y_driver = dat_modis(:); 
         end
         
         Hawaii_dNd_1D_PDF    
         title(['mean=' num2str(Y_mean_overall) ', SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
             ', \DeltaN_d>' num2str(thresh_dNd)]);
        
        
         
     case 'PI and PD'
         switch var_mean
             case {'total_max_random_cloud_amount'}
                 Ybins_DRIVER = [-0.025:0.05:1.025];
                 ylabelstr='Total Cloud Fraction (max random overlap)';
             case {'LWP'}
                 Ybins_DRIVER = [5:10:500];
                 ylabelstr='LWP_{cloudy sky} (g m^{-2})';
             case {'LWPic_using_CF_subgrid_vert_int'}
                 Ybins_DRIVER = [-5:10:500];
                 ylabelstr='LWP_{cloudy sky subgrid} (g m^{-2})';    
             case 'Nd'
                 Ybins_DRIVER = [5:10:500];
                 ylabelstr='Nd (cm^{-3})';
             case 'Nd_using_CF_subgrid_vert_int'
                 Ybins_DRIVER = [5:10:500];
                 ylabelstr='Nd (cm^{-3})';    
             case 'CF'
                 Ybins_DRIVER = [-0.05:0.01:1.05];
                 ylabelstr='CF';
             case 'CF_subgrid_vert_int'
                 Ybins_DRIVER = [-0.05:0.01:1.05];
                 Ybins_DRIVER = [-0.025:0.05:1.025];
                 ylabelstr='CF_{subgrid}';    
             case 'Cloud base rain rates'
                 Ybins_DRIVER = [10.^[-5:0.1:1]];
                 ylabelstr='Cloud base rain rate (mm hr^{-1})';
                 logbin_norm_driver=1; %Normalise by log of bin widths so that the area under the line in log-space corresponds to the 
                 %no. datapoints.
              case 'RWP'
                 Ybins_DRIVER = 10.^[-4:0.1:5];
                 %Ybins_DRIVER = [0:100:;
                 ylabelstr='RWP (g m^{-2})';
                 logbin_norm_driver=1; %Normalise by log of bin widths so that the area under the line in log-space corresponds to the 
                 %no. datapoints.                     
             case 'Liquid Cloud Top Height'
                 %Ybins_DRIVER = 10.^[-4:0.1:5];
                 Ybins_DRIVER = [0:0.2:6];
                 ylabelstr='Liquid Cloud Top Height (km)';
                 %logbin_norm_driver=1; %Normalise by log of bin widths so that the area under the line in log-space corresponds to the 
                 %no. datapoints. 
             case 'Total Ice Water Path'
                 %Ybins_DRIVER = 10.^[-4:0.1:5];
                 %Ybins_DRIVER = 10.^[-8:0.05:0.5];
                 Ybins_DRIVER = 10.^[-7:0.05:0.5]; logbin_norm_driver=1; %Normalise by log of bin widths so that the area under the line in log-space corresponds to the 
                 %no. datapoints. 
                 %Ybins_DRIVER = [-0.025:0.005:1];
                 ylabelstr='Total Ice Water Path (g m^{-2})';                
                 
             case 'Ice Water Path Section 30'
                 %Ybins_DRIVER = 10.^[-4:0.1:5]; logbin_norm_driver=1;
                 Ybins_DRIVER = 10.^[-7:0.05:0.5]; logbin_norm_driver=1; %Normalise by log of bin widths so that the area under the line in log-space corresponds to the 
                 %no. datapoints. 
                 %Ybins_DRIVER = [-0.025:0.005:1];
                 ylabelstr='Ice Water Path Section 30 (g m^{-2})';
                 
                 
         end
         Y_driver = dat_SO2_Nd_PI(:);
         Hawaii_dNd_1D_PDF
         ydat_PI = ydat;
         xdat_PI = xdat;
         mean_from_PDF_PI = Y_mean_overall;
         
         
         Y_driver = dat_SO2_Nd_PD(:);
         Hawaii_dNd_1D_PDF
         ydat_PD = ydat;
         xdat_PD = xdat;
         mean_from_PDF_PD = Y_mean_overall;
         
         figure('color','w'); hold on    
         set(gca,'fontsize',16);
         clear leg_str
         plot(xdat_PI(1).x,ydat_PI(1).y,'b','linewidth',3); leg_str{1}=['Volcano OFF, \mu=' num2str(mean_from_PDF_PI,'%.3g')];
         plot(xdat_PD(1).x,ydat_PD(1).y,'r','linewidth',3); leg_str{2}=['Volcano ON, \mu=' num2str(mean_from_PDF_PD,'%.3g')];         
         title(['SO_2>' num2str(thresh_SO2,'%1.0e') ', LWP>' num2str(thresh_LWP) ...
             ', \DeltaN_d>' num2str(thresh_dNd) ' ' orog_str]); %run_set
         
         legend(leg_str);
         xlabel(ylabelstr);
         if logbin_norm_driver==1
            set(gca,'xscale','log'); 
            ylabel('df/dlog_{10}(x)');
         else
             ylabel('df/dx');
         end
         %set(gca,'xlim',[-5 500]);
         
         
         
 end
 
 
 
 
 
 
 
 %% 2d histogram of LWP volc ON vs OFF
 
X_driver = LWP_PI_ALL; %X_driver(LWP_PI_ALL<10)=NaN;
Y_driver = LWP_PD_ALL; %Y_driver(LWP_PD_ALL<10)=NaN;

Xbins_DRIVER = [0:10:300];
Ybins_DRIVER = Xbins_DRIVER;


    %Y_driver(inan)=NaN;
    %Z_driver = forcing(inds_PI_clear_low);
    xlabelstr='LWP Volcano OFF g m^{-2})';
    ylabelstr = 'LWP Volcano ON (g m^{-2})';
    
    DRIVER_template_2D_histo_LWP_on_vsLWP_off
    shading faceted; %this adds grid lines to the plot - but gives NaNs a
    %colour...
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 100]);
    %label_text_pcolor; %Script to add text labels of the numbers for each block
    title('LWP volc on vs volc off 2D histogram');    

 
 
 %% Cloud base rain rates
 
    vars_UM = {'Rain_rate_cloud_base_in_cloud_LWC_0.05','Rain_rate_cloud_base_in_cloud_LWC_0.01'};          
    
    Hawaii_plot_UM_vars_diffs_etc %Loads the variables and plots PD minus PI diff, PI only and % diff (time
%means).

 
 
 
 
 %% 2d histogram of Rain rate vs LWP/Nd for volc ON

 %LWP is in g/m2, Nd in per m3
 


Nd_thresh = 1e6; %min Nd in per m3
%Nd_thresh = 5e6; %min Nd in per m3
X_driver = LWP_PD_ALL./(Nd_PD_ALL/1e6); X_driver(Nd_PD_ALL<Nd_thresh)=NaN; %Ignore Nd less than 1 per cc
Y_driver = 3600*RainRate_0pt05_PD_ALL; %Y_driver(LWP_PD_ALL<10)=NaN; %convert 

%Xbins_DRIVER = [0:100:4000];
%Ybins_DRIVER = [0:10:350];

Xbins_DRIVER = 10.^[-2:0.02:2];
Ybins_DRIVER = 10.^[-5:0.05:1];

    
    xlabelstr='LWP/Nd (g m^{-2} cm^{3})';
    ylabelstr = 'Rain rate (mm hr^{-1})';
        
    DRIVER_template_2D_histo_RR_vsLWP_over_Nd
    %shading faceted; %this adds grid lines to the plot - but gives NaNs a
    %colour...
    shading flat; %no grid lines
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 100]);
    %label_text_pcolor; %Script to add text labels of the numbers for each block    
    title(['Volcano ON, ' orog_str ', ' emission_str]);
    
    
    rr_Comstock  = 0.0156.*Xbins_DRIVER.^1.75;
    hold on
    plot(Xbins_DRIVER,rr_Comstock,'k');


%% 2d histogram of Rain rate vs LWP/Nd for volc OFF

 %LWP is in g/m2, Nd in per m3
 


Nd_thresh = 1e6; %min Nd in per m3
%Nd_thresh = 5e6; %min Nd in per m3
X_driver = LWP_PI_ALL./(Nd_PD_ALL/1e6); X_driver(Nd_PI_ALL<Nd_thresh)=NaN; %Ignore Nd less than 1 per cc
Y_driver = 3600*Rain_rate_cloud_base_in_cloud_LWC_0pt05; %Y_driver(LWP_PD_ALL<10)=NaN; %convert 


%Bins from Murakami, J. Clim. (2021)
Xbins_DRIVER = 10.^[-2:0.02:2];
Ybins_DRIVER = 10.^[-5:0.05:1];

% Xbins_DRIVER = 10.^[-6:0.02:-1];
% Ybins_DRIVER = 10.^[-5:0.05:1];

    
    xlabelstr='LWP/Nd (g m^{-2} cm^{3})';
    ylabelstr = 'Rain rate (mm hr^{-1})';
        
    DRIVER_template_2D_histo_RR_vsLWP_over_Nd
    %shading faceted; %this adds grid lines to the plot - but gives NaNs a
    %colour...
    shading flat; %no grid lines
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([0 100]);
    %label_text_pcolor; %Script to add text labels of the numbers for each block
    title(['Volcano OFF, ' orog_str ', ' emission_str]);
    
    
    rr_Comstock  = 0.0156.*Xbins_DRIVER.^1.75;
    hold on
    plot(Xbins_DRIVER,rr_Comstock,'k');
    
    
    
%% Experiments with LWP threshold for cloud fraction. 
thresh_LWP = [0:5:150];
thresh_LWP = [0:1:50];
thresh_LWP = [0:0.01:1];

%Also filter using SO2 to restrict to plume region

%SO2 thresholds
thresh_SO2 = -1;
%thresh_SO2 = 1e-5;
%thresh_SO2 = 1e-4;

var_test_01 = LWP_PI_ALL;
var_test_02 = LWP_PD_ALL;

inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);        
var_test_01(inan)=NaN; 
var_test_02(inan)=NaN; 

clear N_PI N_PD leg_str
for i=1:length(thresh_LWP)
    N_PI(i) = length(find(var_test_01 >= thresh_LWP(i)));
    N_PD(i) = length(find(var_test_02 >= thresh_LWP(i)));    
end

Ntot = length(LWP_PI_ALL(:)) - length(inan);
CF_PI = N_PI/Ntot; leg_str{1} = 'Volcano OFF';
CF_PD = N_PD/Ntot; leg_str{2} = 'Volcano ON';

figure('color','w');
plot(thresh_LWP,CF_PI,'linewidth',3); hold on
plot(thresh_LWP,CF_PD,'r','linewidth',3);
xlabel('LWP threshold (g m^{-2})');
ylabel('Fraction of datapoints');
legend(leg_str);
fontsize_figure(gcf,gca,18); 

%% As above, but using cloud optical depth

%Estimate COD using LWP and Nd.

CTT=278; P=850e2; k=0.8; 
[tau_PI_ALL,H_PI]=MODIS_tau_func_N_LWP(LWP_PI_ALL*1e-3,Nd_PI_ALL,CTT,k,P); %LWP needs to be in kg/m2 and Nd in per m3
[tau_PD_ALL,H_PD]=MODIS_tau_func_N_LWP(LWP_PD_ALL*1e-3,Nd_PD_ALL,CTT,k,P); 

dx=0.1; thresh_tau = [0:dx:10];
dx=0.01; thresh_tau = [0:dx:1]; 

%SO2 thresholds
thresh_SO2 = -1;
thresh_SO2 = 1e-5;
%thresh_SO2 = 1e-4;

var_test_PI = tau_PI_ALL;
var_test_PD = tau_PD_ALL;

inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);        
var_test_PI(inan)=NaN; 
var_test_PD(inan)=NaN; 

clear N_PI2 N_PD2 leg_str
for i=1:length(thresh_tau)
    N_PI2(i) = length(find(var_test_PI >= thresh_tau(i)));
    N_PD2(i) = length(find(var_test_PD >= thresh_tau(i)));    
end

Ntot_type = 'ALL';
%Ntot_type = 'Individual';

switch Ntot_type
    case 'ALL'
        Ntot_PI = length(tau_PI_ALL(:)) - length(inan); %total number minus those outside of the plume area
        Ntot_PD = Ntot_PI;
    case 'Individual'
        %Ntot = length(tau_PI_ALL(:)) - length(inan); %total number minus those outside of the plume area
        Ntot_PI = length(find(var_test_PI > -1));
        Ntot_PD = length(find(var_test_PD > -1));
end


%(since we only only considering the plume region).
CF_PI2 = N_PI2/Ntot_PI; leg_str{1} = 'Volcano OFF';
CF_PD2 = N_PD2/Ntot_PD; leg_str{2} = 'Volcano ON';

figure('color','w');
plot([thresh_tau],CF_PI2,'linewidth',3); hold on
plot([thresh_tau],CF_PD2,'r','linewidth',3);
xlabel('\tau threshold');
ylabel('Fraction of datapoints');
legend(leg_str);
title(['SO_2 threshold = ' num2str(thresh_SO2)]);
fontsize_figure(gcf,gca,18); 


%% Test using the same population of clouds, but with Nd increased by a certain factor
%Estimate COD using LWP and Nd.

CTT=278; P=850e2; k=0.8; 
[tau_PI_ALL,H_PI]=MODIS_tau_func_N_LWP(LWP_PI_ALL*1e-3,Nd_PI_ALL,CTT,k,P); %LWP needs to be in kg/m2 and Nd in per m3
fNd = 2; %Double Nd
[tau_PD_ALL,H_PD]=MODIS_tau_func_N_LWP(LWP_PI_ALL*1e-3,Nd_PI_ALL*fNd,CTT,k,P); %LWP needs to be in kg/m2 and Nd in per m3


%dx=0.1; thresh_tau = [-dx dx:dx:10];
dx=0.1; thresh_tau = [0:dx:10];
dx=0.01; thresh_tau = [0:dx:1]; 

%SO2 thresholds
%thresh_SO2 = -1;
thresh_SO2 = 1e-5;
%thresh_SO2 = 1e-4;



var_test_PI = tau_PI_ALL;
var_test_PD = tau_PD_ALL;

inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);        
var_test_PI(inan)=NaN; 
var_test_PD(inan)=NaN; 

clear N_PI2 N_PD2 leg_str
for i=1:length(thresh_tau)
    N_PI2(i) = length(find(var_test_PI >= thresh_tau(i)));
    N_PD2(i) = length(find(var_test_PD >= thresh_tau(i)));    
end


Ntot_type = 'ALL';
Ntot_type = 'Individual';

switch Ntot_type
    case 'ALL'
        Ntot_PI = length(tau_PI_ALL(:)) - length(inan); %total number minus those outside of the plume area
        Ntot_PD = Ntot_PI;
    case 'Individual'
        %Ntot = length(tau_PI_ALL(:)) - length(inan); %total number minus those outside of the plume area
        Ntot_PI = length(find(var_test_PI > -1));
        Ntot_PD = length(find(var_test_PD > -1));
end

%(since we only only considering the plume region).
CF_PI2 = N_PI2/Ntot_PI; leg_str{1} = 'Volcano OFF';
CF_PD2 = N_PD2/Ntot_PD; leg_str{2} = ['Volcano OFF, Nd*' num2str(fNd)];

figure('color','w');
plot([thresh_tau],CF_PI2,'linewidth',3); hold on
plot([thresh_tau],CF_PD2,'r','linewidth',3);
xlabel('\tau threshold');
ylabel('Fraction of datapoints');
legend(leg_str);
title(['SO_2 threshold = ' num2str(thresh_SO2)]);
fontsize_figure(gcf,gca,18); 



%% Test using the same population of clouds, but with Nd increased by a certain factor - CF vs Nd factor for chosen tau thresh
%Estimate COD using LWP and Nd.

CTT=278; P=850e2; k=0.8; 
[tau_PI_ALL,H_PI]=MODIS_tau_func_N_LWP(LWP_PI_ALL*1e-3,Nd_PI_ALL,CTT,k,P); %LWP needs to be in kg/m2 and Nd in per m3
fNds = [0.2:0.2:5]; %Double Nd
clear tau_fNd
for i=1:length(fNds)
    [tau_fNd{i},H_PD]=MODIS_tau_func_N_LWP(LWP_PI_ALL*1e-3,Nd_PI_ALL*fNds(i),CTT,k,P); %LWP needs to be in kg/m2 and Nd in per m3
end

%dx=0.1; thresh_tau = [-dx dx:dx:10];
%dx=0.1; thresh_tau = [0:dx:10];
%dx=0.01; thresh_tau = [0:dx:1];

thresh_tau = 0.3; %The MODIS threshold I think
thresh_tau = 10; %

%SO2 thresholds
%thresh_SO2 = -1;
thresh_SO2 = 1e-5;
%thresh_SO2 = 1e-4;

inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);        

clear N_PI2 N_PD2 leg_str
for i=1:length(tau_fNd)
    
    var_test_PI = tau_fNd{i};
    %var_test_PD = tau_PD_ALL;    
    
    var_test_PI(inan)=NaN;
    %var_test_PD(inan)=NaN;    
    
    N_PI2(i) = length(find(var_test_PI >= thresh_tau));
    %N_PD2(i) = length(find(var_test_PD >= thresh_tau(i)));
end


%Ntot = length(tau_PI_ALL(:)) - length(inan); %total number minus those outside of the plume area

Ntot_PI = length(find(var_test_PI > -1));
%Ntot_PD = length(find(var_test_PD > -1));

CF_PI2 = N_PI2/Ntot_PI; leg_str{1} = ['\tau threshold=' num2str(thresh_tau)];
%CF_PD2 = N_PD2/Ntot_PD; leg_str{2} = ['Volcano OFF, Nd*' num2str(fNd)];

figure('color','w');
plot(fNds,CF_PI2,'linewidth',3); hold on
%plot([thresh_tau],CF_PD2,'r','linewidth',3);
xlabel('Nd factor');
ylabel('Fraction of datapoints');
legend(leg_str);
title(['SO_2 threshold = ' num2str(thresh_SO2)]);
fontsize_figure(gcf,gca,18); 
grid on




