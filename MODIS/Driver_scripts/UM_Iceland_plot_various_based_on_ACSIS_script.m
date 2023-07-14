%Pre-process using
%   UM_quick_plot_global.m

%Add variables to climits are set here too :-
%   UM_var_defs.m

%runs this plotting script - lat/lon for map also set here. FOR PDFs, etc. it is set below :-
  %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
  
% Here, PI will refer to volcano OFF, PD volcano ON

UM_ACSIS_SW_vs_cloud_properties_global_DEFAULTS

savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/Iceland/plots_' datestr(now,30) '/'];
eval(['!mkdir ' savedir_date]);

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

%% Choose runs
run_set = 'global runs Sep 2018';
run_set = 'nested UKCA runs Sep 2018'; %u-ba333 and u-ba050
%run_set = 'nested UKCA-CASIM runs Sep 2018';
%run_set = 'nested UKCA-CASIM-fixed_sed runs Nov 2018'; %u-bc309 and u-bc308
%run_set = 'nested UKCA vs CASIM, volc OFF runs Sep 2018';
%run_set = 'nested UKCA vs CASIM, volc ON runs Sep 2018';

Nd_var_str='UKCA';

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

ioverride_LAT_plots=1;
LAT_val_DRIVER_override = [40 82]; LON_val_DRIVER_override = [-85 55]; %N Atlantic Iceland region


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
nest=load(filename_nest);
irotated_pole_box=0;


%Choose a variable to use to get the times
%var_UM_DRIVER = 'LWP';
%var_UM_DRIVER = 'accum_number_ukca';
%var_UM_DRIVER = 'SW_down_surf';
%var_UM_DRIVER = 'Nd_lwc_weighted_UKCA';
var_UM_DRIVER = 'LWP_time_only'; %use the LWP file, but only for loading the time
var_UM_DRIVER = 'LS_surf_rain_rate_time_only'; %use the file, but only for loading the time

scrsz=get(0,'ScreenSize');
posit=[9 60 scrsz(3) scrsz(4)];
figure('position',posit);

%% Choose times/dates here - doing all times for now
%date_str_range = 'all';

% time_round = time_ALL(it_global_diff); %this gets put on the plots
% time_round = datenum(date_str); %this gets put on the plots

% dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
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

%Specific days :-
%clear time_choice; time_choice.time_specific = datenum('06-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('07-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('10-Aug-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Dec-2009'); time_choice.find_nearest=1;





%% Get time data from the specified file
    var_UM = var_UM_DRIVER;
    
    switch run_type_DRIVER
        case 'nested'
            um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'nested'; icoarse=0; ivar_dir=1; %nested UKCA
            dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];

        case 'global'
            % PI run
            %    um_case='u-au652'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
            um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
            dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];

    end

    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

    

array_in=[]; %just test to get the indices for now.
dim=NaN; %don't need dim if just getting the indices
[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);
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
    %[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
    [gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);    
    %    gcm_Plon2D_edges_UM = dat_global.gcm_Plon2D_edges_UM;

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
 
um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

array_in=[]; %just test to get the indices for now.
dim=NaN; %don't need dim if just getting the indices
[out, time_out_SW_in, time_inds_SW_in, dtime_match_SW_in] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

[SW_in_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_SW_in,load_type,gcm_Plat2D_UM);

um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
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
 
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [SW_down_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
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
 
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [SW_nodirect_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [SW_nodirect_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);    
    
    aerosol_effect_PI_ALL = SW_down_PI_ALL - SW_nodirect_PI_ALL;    
    aerosol_effect_PD_ALL = SW_down_PD_ALL - SW_nodirect_PD_ALL;
    direct_ALL = aerosol_effect_PD_ALL - aerosol_effect_PI_ALL;
    
    
% For indirect effect we do clean minus clear_clean.
% I.e., no aerosol effect with cloud effect - no aerosol effect without cloud effect
    %no direct or indirect effects
    
 var_UM = [SW_var '_clean_clear_'  SW_surf_or_TOA]; 
    
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [SW_clean_clear_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
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

%% Change in SW out TOA and SW out TOA vs CERES
isave_plot_CERES=1;

    clear gca
    var_UM = 'SW_up_TOA';
 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [SW_TOA_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [SW_TOA_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    dat_modis = meanNoNan(SW_TOA_PD_ALL,3) - meanNoNan(SW_TOA_PI_ALL,3); var_UM = 'Change in SW out TOA (PD minus PI)';
    subtitle_str = var_UM;
%    dat_modis = meanNoNan(SW_TOA_PD_ALL,3);
%    dat_modis = meanNoNan(SW_TOA_PI_ALL,3);
   
    
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
 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [low_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    low_CF_PD_ALL(low_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??  
    
    if icosp_mask==1        
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [low_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        low_CF_PD_ALL(low_CF_mask==0) = NaN;
    end
    
    var_UM = 'low_cloud_amount';
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [low_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    low_CF_PI_ALL(low_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??
    
    low_CF_av_ALL = (low_CF_PI_ALL + low_CF_PD_ALL) /2;
        
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

    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [mid_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    mid_CF_PD_ALL(mid_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [mid_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        mid_CF_PD_ALL(mid_CF_mask==0) = NaN;
    end
    
%Don't have COSP output for PI yet
    var_UM = 'mid_cloud_amount';
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);     
    [mid_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    mid_CF_PI_ALL(mid_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??   
    
    mid_CF_av_ALL = (mid_CF_PI_ALL + mid_CF_PD_ALL) /2;
    

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
    
 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [high_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    high_CF_PD_ALL(high_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [high_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        high_CF_PD_ALL(high_CF_mask==0) = NaN;
    end
    
    var_UM = 'high_cloud_amount';
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type); 
    [high_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    high_CF_PI_ALL(high_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??   
    
    high_CF_av_ALL = (high_CF_PI_ALL + high_CF_PD_ALL) /2;
       
    
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

        um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        [total_CF_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??

        if icosp_mask==1
            dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
            [total_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
            total_CF_PD_ALL(total_CF_mask==0) = NaN;
        end

        %     var_UM = 'total_cloud_amount';
        %     um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        %     dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        %     dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
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
iload_modis_dat=1;    

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
 %Nd_type = 'lwc in-cloud weighted'; 
 
    switch Nd_type
        case 'lwc weighted'
            var_UM = ['Nd_lwc_weighted_' Nd_var_str];
        case 'lwc in-cloud weighted'
            var_UM = ['Nd_lwc_in_cloud_weighted_' Nd_var_str];        
    end
    
    icosp_mask=0;
    
        
 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [Nd_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
%    total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [total_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        total_CF_PD_ALL(total_CF_mask==0) = NaN;
    end
    

    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [Nd_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
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
    
%Change in Nd from PI to PD

    Nd_PI_map = meanNoNan(Nd_PI_ALL,3);
    Nd_PD_map = meanNoNan(Nd_PD_ALL,3);
    
    dat_modis = (Nd_PD_map - Nd_PI_map)/1e6; var_UM = 'Change in Nd (PD minus PI; cm^{-3})';    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]);
    
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
    
%% Plot LWP
isave_plot_global_LWP=0;
isave_plot_LWP_bias=0;
iplot_amsre=0;

clear gca

 var_UM_LWP = 'LWP'; %calculated from Python script using QCL after timestep(0-254)
% var_UM = 'LWP_sec30'; %from section 30-405

    var_UM = var_UM_LWP;
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [LWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    
    var_UM = var_UM_LWP; 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [LWP_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);   
    
    
    switch var_UM
        case {'LWP_sec30','LWP_correct_dz'}
           LWP_PD_ALL =  LWP_PD_ALL*1e3; %convert to g/m2
           LWP_PI_ALL =  LWP_PI_ALL*1e3; %convert to g/m2           
    end
    
    mean_LWP_PI = meanNoNan(LWP_PI_ALL,3);
    mean_LWP_PD = meanNoNan(LWP_PD_ALL,3);  
    
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
    set(gca,'ylim',[-2 8])
    fontsize_figure(gcf,gca,18);
    grid on
    
   
    
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
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [LS_rain_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    
    var_UM = var_UM_LWP; 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
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

dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case_PD];
save_dLWP_vs_rain_rates = [dirUM '/save_dLWP_vs_rain_rates.mat'];

save(save_dLWP_vs_rain_rates,'-V7.3','Y_me_time_mean','Y_N_time_mean','Y_std_time_mean','LS_rain_bins_time_mean');
save(save_dLWP_vs_rain_rates,'-APPEND','Y_me','Y_N','Y_std','LS_rain_bins','mean_LS_rain_PI','mean_LS_rain_PD');





%% RWP, convective LWP, etc.
    
    iload_rwp=0;
    if iload_rwp==1
        var_UM = 'RWP';
        um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        [RWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);


        var_UM = 'RWP';
        um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
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
        um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
   
        %The fisrt time is missing for the convective LWP diags (5-213 etc)
        %for some reason - so try using a special time indices for this, since the time fields are different.

        array_in=[]; %just test to get the indices for now.
        dim=NaN; %don't need dim if just getting the indices
        [out, time_out_conv_LWP, time_inds_conv_LWP, dtime_match_conv_LWP] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);
        
        [Conv_LWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_conv_LWP,load_type,gcm_Plat2D_UM);


        var_UM = 'Conv_LWP';
        um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        [Conv_LWP_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_conv_LWP,load_type,gcm_Plat2D_UM);

        dat_modis = meanNoNan(Conv_LWP_PD_ALL,3);

        %run plotting script
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        caxis([0 150]);
        
        
        var_UM = 'Conv_RWP';
        um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);          
        [Conv_RWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_conv_LWP,load_type,gcm_Plat2D_UM);


        var_UM = 'Conv_RWP';
        um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
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
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [RR_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    
    var_UM = var_UM_RR; 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
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
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    
    var_UM = 'low_cloud_amount';
    
%     filename = [dirUM '/' var_UM '/' run_type '_' var_UM '_native res_ALL.mat'];
%     dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

    
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


