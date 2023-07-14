% N.B. - to get a global plot can just set irestrict_domain_DRIVER=0

%Pre-process using
%   UM_quick_plot_global.m

%Add variables to climits are set here too :-
%   UM_var_defs.m

%runs this plotting script - lat/lon for map also set here. FOR PDFs, etc. it is set below :-
  %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
  
UM_ACSIS_SW_vs_cloud_properties_global_DEFAULTS
  
savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_' datestr(now,30) '/'];
%savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_TEMP/'];
eval(['!mkdir ' savedir_date]);  

cf_thresh=0.01;

i_Iceland=0;
icoarse_grain=0;
M_coarse_grain=4; N_coarse_grain=4;
  
close_all_figs=1;  
  
load_type = 'mat';
load_type = 'merged netCDF';

cloud_input = 'UM'; %Use the usual UM low, mid and high cloud fractions.
%cloud_input = 'CALIPSO'; %Use the UM COSP CALIPSO values
%cloud_input = 'MODIS'; %Use the COSP MODIS values


var_UM_DRIVER = 'LWP';
var_UM_DRIVER = 'accum_number_ukca';
var_UM_DRIVER = 'SW_down_surf';
%var_UM_DRIVER = 'Nd_lwc_weighted_UKCA';
var_UM_DRIVER = 'LWP_sec30';

%[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

iload_rwp=1;
iload_SW_down_TOA=1;
amsre_data = 1;

ical_data=1; %Whether are using CALIPSO data - load using script :- read_calipso_monthly_night_IPSL.m for average values
    %first. This puts the data into these fields :- 
    % cllcalipso_monthly_AVERAGE, clmcalipso_monthly_AVERAGE, clhcalipso_monthly_AVERAGE
    % Is average of day and night values for now (have separate ones too).

iceres_data=1; %whether to load CERES data
icontour_DRIVER=0;

%run_set = 'sst updating, temperature nudging';
run_set = 'wind only nudging from level 18';
%run_set = 'NO AIEs, wind only nudging from level 18';
%run_set = 'NO rain AIEs, wind only nudging from level 18'
run_set = 'wind only nudging from level 18, reverted Hamish switches';
run_set = 'wind only nudging from level 18, reverted Hamish switches, aerosol sw bug fix';
%run_set = 'wind only nudging from level 18, reverted Hamish switches, aerosol sw bug fix, KOGAN';
%run_set = 'CAM6 first set of output 4th Feb 2019';
%run_set = 'wind only nudging from level 18, reverted Hamish switches, aerosol sw bug fix, 7 hour output RyanE';
%run_set = 'Rosenfeld, no AIEs'; %(Rosenfeld runs) - workspace - save_file2='/home/disk/eos10/d.grosvenor/UM/Rosenfeld_analysis/Rosenfeld_workspace.mat';
%run_set = 'Rosenfeld, no AIEs, no aerosol scavenging'; %Only ran to start of 14th May, 2010 before BicStab
%run_set = 'Rosenfeld, no AIEs, no aerosol scavenging, NUDGING';        
%run_set = 'Rosenfeld, no AIEs, scavenging ON, NUDGING';
%run_set = 'Rosenfeld, no AIEs, no convective scavenging vs all scavenging with bug fix, NUDGING';
%run_set = 'Rosenfeld, no AIEs, scavenging ON with bug fix, NUDGING';
%run_set = 'Rosenfeld, no AIEs, scavenging ON inc. convective with bug fix, NUDGING';
run_set = 'Ken/Leighton tests of scavenging params, NUDGING';

switch run_set
    case 'no sst updating'
        um_case_PI = 'u-av504';        
        um_case_PD = 'u-av503'; %No SST updating

    case 'sst updating, temperature nudging'
        um_case_PI = 'u-ay364';
        um_case_PD = 'u-ax981';

    case 'wind only nudging from level 18'
        um_case_PI = 'u-ay837';
        um_case_PD = 'u-ay838';
        
    case 'NO AIEs, wind only nudging from level 18'
        um_case_PI = 'u-bc437';
        um_case_PD = 'u-bc438';
        
    case 'NO rain AIEs, wind only nudging from level 18'
        um_case_PI = 'u-bc675';
        um_case_PD = 'u-bc676';        
        
    case 'wind only nudging from level 18, reverted Hamish switches'
        um_case_PI = 'u-be716';
        um_case_PD = 'u-be717';
        
    case 'wind only nudging from level 18, reverted Hamish switches, aerosol sw bug fix'
        um_case_PI = 'u-bf109';
        um_case_PD = 'u-bf110';
        
     case 'wind only nudging from level 18, reverted Hamish switches, aerosol sw bug fix, KOGAN'
        um_case_PI = 'u-bf169';
        um_case_PD = 'u-bf170';
        
    case 'CAM6 first set of output 4th Feb 2019'
        um_case_PI = 'cam6_PI_04Feb2019';
        um_case_PD = 'cam6_PD_04Feb2019';
        load_type = 'cam6';
        var_UM_DRIVER='FSDS';
        iload_SW_down_TOA=0;
        ical_data=0;
        iload_rwp=0;
        
    case 'Rosenfeld, no AIEs'
        um_case_PI = 'u-bf666';
        um_case_PD = 'u-bf666';
        
    case 'wind only nudging from level 18, reverted Hamish switches, aerosol sw bug fix, 7 hour output RyanE'
        um_case_PI = 'u-bi042';
        um_case_PD = 'u-bi042';
        
    case 'Rosenfeld, no AIEs, no aerosol scavenging'
        um_case_PI = 'u-bf666';
        um_case_PD = 'u-bi194';    
        
    case 'Rosenfeld, no AIEs, no aerosol scavenging, NUDGING'
        um_case_PI = 'u-bi724';
        um_case_PD = 'u-bi725';
        
    case 'Rosenfeld, no AIEs, scavenging ON, NUDGING'
        um_case_PI = 'u-bi725';
        um_case_PD = 'u-bi724';
                
    case 'Rosenfeld, no AIEs, no convective scavenging vs all EXCEPT CONVECTIVE scavenging with bug fix, NUDGING'
        %first year of run :-    
        %'/home/disk/eos10/d.grosvenor/UM/Rosenfeld_analysis/Rosenfeld_scav_OFF_plus_noconv_workspace_nudging.mat'
        um_case_PD = 'u-bi989'; %Scavenging OFF (including convective/plume scavenging)
        um_case_PI = 'u-bi990'; %Scavenging ON (scavenging bug fixed)
        um_case_CONV_SCAV_ON = 'u-bi945'; %Convective scavenging ON, other scavengign OFF.

        
    case 'Rosenfeld, no AIEs, scavenging ON with bug fix, NUDGING'
    %first year of run :-    
      %/home/disk/eos10/d.grosvenor/UM/Rosenfeld_analysis/Rosenfeld_scav_ON_bugfix_workspace_nudging.mat
    %last year of run (28-Mar-2011 to 28-Mar-2012) :-
      %/home/disk/eos10/d.grosvenor/UM/Rosenfeld_analysis/Rosenfeld_scav_ON_bugfix_workspace_nudging_last_year.mat        
        um_case_PD = 'u-bi990'; %Scavenging ON (scavenging bug fixed)
        um_case_PI = 'u-bi990'; %Scavenging ON (scavenging bug fixed)
        um_case_CONV_SCAV_ON = 'u-bi945'; %Convective scavenging ON, other scavengign OFF.   
        
  case 'Rosenfeld, no AIEs, scavenging ON inc. convective with bug fix, NUDGING'
        %first year of run :-    
        %'/home/disk/eos10/d.grosvenor/UM/Rosenfeld_analysis/Rosenfeld_scav_ON_inc_conv_u_bf555_workspace_nudging.mat'        
        um_case_PI = 'u-bi989'; %Scavenging OFF (including convective/plume scavenging)
        um_case_PD = 'u-bj555'; %Scavenging ON (scavenging bug fixed, inc. convective)
        um_case_CONV_SCAV_ON = 'u-bi945'; %Convective scavenging ON, other scavenging OFF.
        um_case_CONV_SCAV_OFF = 'u-bi990'; %Convective scavenging OFF, other scavenging ON. - not using this yet
        % So, should do ALL scav (u-bj555) vs conv only OFF (u-bi990) and
        % ALL scav vs other scav OFF (u-bi945)
        
  case 'Ken/Leighton tests of scavenging params, NUDGING'
        um_case_PI = 'u-bi989'; %Scavenging OFF (including convective/plume scavenging)
        um_case_PD = 'u-bj555'; %Scavenging ON (scavenging bug fixed, inc. convective)
        um_case_CONV_SCAV_ON = 'u-bi945'; %Convective scavenging ON, other scavenging OFF.
        um_case_CONV_SCAV_OFF = 'u-bi990'; %Convective scavenging OFF, other scavenging ON. - not using this yet
        
        
end

time_format_str=' UTC';

%% set the box region within which to do the cloud state analysis, timeseries, etc.
%% N.B. the lat/lon for maps is set in the following script (not by LAT_val_DRIVER2 below) :-
  %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
  % N.B. - to get a global plot can just set irestrict_domain_DRIVER=0


region_choice = 'Southern NA';
%region_choice = 'Northern NA';
region_choice = 'Rosenfeld VOCALS';
region_choice = 'Rosenfeld ALL';
%region_choice = 'VOCALS CPT';
%region_choice = 'VOCALS coastal';

% Run the script to pick the region
[LAT_val_DRIVER2, LON_val_DRIVER2, region_shortname] = UM_ACSIS_choose_region(region_choice);


%% Other switches etc.

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

irestrict_domain_DRIVER=1;  %N.B. - have to re-run some of the commands below when using this...

clear gca

%If want to plot an outline of the nest on the global map
filename_nest = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/BC__model_level01_time08.mat';
nest=load(filename_nest);
irotated_pole_box=0;




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

clear time_choice; time_choice.time_range = [datenum('28-Mar-2009') datenum('29-Mar-2010')];
%clear time_choice; time_choice.time_range = [datenum('01-Apr-2009') datenum('31-Mar-2012')]; %full time period of Rosenfeld runs
%clear time_choice; time_choice.time_range = [datenum('28-Mar-2011') datenum('29-Mar-2012')]; %Rosenfeld run u-bf666
%clear time_choice; time_choice.time_range = [datenum('01-Mar-2009') datenum('29-Mar-2010')];
%clear time_choice; time_choice.time_range = [datenum('28-Mar-2009') datenum('01-Nov-2009')];
%clear time_choice; time_choice.time_range = [datenum('01-Aug-2009') datenum('31-Aug-2009')];
%clear time_choice; time_choice.time_range = [datenum('01-Dec-2009') datenum('31-Dec-2009')];

%Specific days :-
%clear time_choice; time_choice.time_specific = datenum('06-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('07-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('10-Aug-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Dec-2009'); time_choice.find_nearest=1;



%% load the driver variable to get time indices etc.
isave_plot_global_diff=0;

scrsz=get(0,'ScreenSize');
posit=[9 60 scrsz(3) scrsz(4)];
figure('position',posit);





if iplot_global==1
    var_UM = var_UM_DRIVER;

% PI run    
%    um_case='u-au652'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];

    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

    

array_in=[]; %just test to get the indices for now.
dim=NaN; %don't need dim if just getting the indices
[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);
time_round = datestr(dat_global.time_ALL(time_inds(1)));

[Y_UM,M_UM,D_UM,HH_UM,MM_UM,SS_UM] = datevec(time_out);


    
    gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
    gcm_Plat2D_edges_UM = dat_global.gcm_Plat2D_edges_UM;

    gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
    
    %if irestrict_domain_DRIVER==0  %Put this code in plot_global_maps now.
    %    i180 = find(gcm_Plon2D_UM<0);
    %    gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180)+360;
    %end
    [gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
    
    
    [temp_var,nT_main] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    nT=nT_main;
    
    
%% global run SW surface forcing map
    
    var_UM = 'SW_down_surf';
    
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
    [SW_down_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);    
    
    dat_modis_PI = meanNoNan(SW_down_PI_ALL,3);
    
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
       
    [iregion_lin,iregion_lin_edges,SW_down_NA_PI]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,SW_down_PI_ALL(:,:,1:end)); 
    SW_down_NA_PI_timser = meanNoNan(meanNoNan(SW_down_NA_PI,1),1);

    
% PD run
%    um_case='u-au536'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    %filename = [dirUM '/' var_UM '/' run_type '_' var_UM '_native res_ALL.mat'];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

    
    array_in=[]; %just test to get the indices for now.
    dim=NaN; %don't need dim if just gettin the indices
    [out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

%    nT = length(eval(['dat_global.' var_UM '_ALL;']));
%     nT = length(time_inds);
%     it=0;
%     SW_down_PD_ALL = NaN * ones([size(gcm_Plat2D_UM,1) size(gcm_Plat2D_UM,2) nT]);    
%     for it_global_diff=time_inds %1:nT
%         it=it+1;
%         SW_down_PD_ALL(:,:,it) = eval(['dat_global.' var_UM '_ALL{it_global_diff};']);        
%     end
    
    [SW_down_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    dat_modis_PD = meanNoNan(SW_down_PD_ALL,3);
    
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
       
    [iregion_lin,iregion_lin_edges,SW_down_NA_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,SW_down_PD_ALL(:,:,1:end)); 
    SW_down_NA_PD_timser = meanNoNan(meanNoNan(SW_down_NA_PD,1),1);
 
    
   



    dat_modis = dat_modis_PD - dat_modis_PI;
    
    %run plotting script
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global



end



% Save
%savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/model_vs_MODIS_Nd_plots/' modis_day_of_year_str '_global_' titlenam_driver];
savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];

if isave_plot_global_diff==1
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
    close(gcf);
end



if iload_SW_down_TOA==1
    SW_down_type = 'SW_down_clean_clear_TOA';
    
    switch SW_down_type
        case 'SW_down_clean_clear_TOA'
            %TOA SW downwards (incoming) - The clean and clean_clear TOA
            %diags are on the same timesteps as the surface ones and so are
            %the better choice
            var_UM = SW_down_type;
            
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
            
        case 'SW_down_TOA'
            
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
    
    
end

%% Quick plot of aerosol number for Hamish
iaerol_Hamish_plots=0;
if iaerol_Hamish_plots==1
    var_UM = 'nucleation_number_ukca';
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
    N_nuc=dat_global.dat;
    
    var_UM = 'aitken_number_ukca';
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
    N_aitken=dat_global.dat;
    
    var_UM = 'accum_number_ukca';
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
    N_accum=dat_global.dat;
    
    var_UM = 'coarse_number_ukca';
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
    N_coarse=dat_global.dat;   
        
    Ntot_all = N_nuc + N_aitken + N_accum + N_coarse;
    Ntot_no_nuc = N_aitken + N_accum + N_coarse;
    
    dat_modis = meanNoNan(Ntot_all,1);
    %run plotting script
    figure
    var_UM = 'Nucleation + Aitken + accumulation + coarse aerosol number on level 39 (10 km, # per # air molecules)'
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 5e-15]);

    dat_modis = meanNoNan(Ntot_no_nuc,1);
    %run plotting script
    figure
    var_UM = 'Aitken + accumulation + coarse aerosol number on level 39 (10 km, # per # air molecules)'
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 5e-16]);
    
end

%% Load satellite data, etc.

 % NOTE - might have to be a little careful with the MODIS data since 
 % each daily could be an average of more than one overpass. The screening
 % for CF is done on each overpass individually. So, data for CF>80 and
 % CF>0 on the same day might not be the same since there could be one
 % overpass that has a CF<80 and so is rejected and one that is kept. For
 % the CF>0 screening they will both be kept and so we can get different
 % CFs, etc.
 % This might affect the co-location of the MODIS and AMSRE data since
 % AMSRE and AMSR2 have swath widths of only 1450km compared to 2300km for
 % MODIS. So, MODIS may average in more than one swath, but not AMSR. Can
 % perhaps guard against this by restricting the sensor angle to only
 % include the AMSR swath width.
    
 % saved some global MODIS (mockL3) data for the period [datenum('28-Mar-2009') datenum('29-Mar-2010')]; here :-
 %CF>80
 modis_Nd_file = '/home/disk/eos10/d.grosvenor/saved_MODIS_data/NAtlantic/global_mockL3_CF80_SZA65_2009-2010_ACSIS.mat';
 %CF>0
 modis_Nd_file_CF0 = '/home/disk/eos10/d.grosvenor/saved_MODIS_data/NAtlantic/global_mockL3_CF0_SZA65_2009-2010_ACSIS.mat'; 
 modis_Nd_file_CF0_2 = '/home/disk/eos10/d.grosvenor/saved_MODIS_data/NAtlantic/global_mockL3_CF0_SZA65_no_CTH_CTT_screening_2009-2010_ACSIS.mat';  
 % Load on first use :-
 if ~exist('Cloud_Fraction_Liquid')
     load(modis_Nd_file,'CTH','W_time3','Cloud_Fraction_Liquid','Droplet_Number_Concentration_37',...         
         'gcm_Plat2D_AMSRE','gcm_Plon2D_AMSRE','modisyear_timeseries3','daynum_timeseries3');          
     
     %Somehow daynum_timeseries3 is getting screwed up... fix this
     %CF>80 Nd data
     modis_loaded =  load(modis_Nd_file,'CTH','W_time3','Cloud_Fraction_Liquid','Droplet_Number_Concentration_37',...
         'gcm_Plat2D_AMSRE','gcm_Plon2D_AMSRE','modisyear_timeseries3','daynum_timeseries3');  
     
     clear opts
     opts.iscreen_seaice = 1;
     opts.iscreen_land = 0; %Ideally need this to remove the ice covered surfaces that may not be sea-ice.
        %But may want to keep land in some cases - I guess would not want
        %to look at land in Antarctica, so could remove it all there, but
        %not in other places?
     opts.iseaice_max_2week=1; 
     opts.mask_size=3; %How much to smooth the land mask - using 3 here gets rid of the ice shelf region near the Antarctic Peninsula
     %opts.vars = {'N_time3','N_time3_37','Droplet_Number_Concentration.timeseries3','Droplet_Number_Concentration_37.timeseries3','Cloud_Optical_Thickness_Liquid_Mean.timeseries3',...
        %'Cloud_Effective_Radius_Liquid_Mean.timeseries3','Cloud_Effective_Radius_37_Liquid_Mean.timeseries3'};
     opts.vars = {'Droplet_Number_Concentration_37.timeseries3','W_time3','Cloud_Fraction_Liquid.timeseries3','CTH.timeseries3'};
     opts.direcs={'aqua'};  opts.years_multi=unique(modis_loaded.modisyear_timeseries3);
     
     
     [modis_loaded]=screen_seaice_MODIS(modis_loaded,opts);
     
     
     
     modis_loaded_CF0 =  load(modis_Nd_file_CF0,'CTH','W_time3','Cloud_Fraction_Liquid','Droplet_Number_Concentration_37',...
         'Cloud_Fraction_NoOpt','Cloud_Fraction_Combined','Cloud_Fraction_Undetermined',...
         'gcm_Plat2D_AMSRE','gcm_Plon2D_AMSRE','modisyear_timeseries3','daynum_timeseries3');      
     %modis_loaded_CF0 =  load(modis_Nd_file_CF0,'CTH','W_time3','Droplet_Number_Concentration_37',...
     %    'gcm_Plat2D_AMSRE','gcm_Plon2D_AMSRE','modisyear_timeseries3','daynum_timeseries3');        
     opts.years_multi=unique(modis_loaded_CF0.modisyear_timeseries3);
     [modis_loaded_CF0]=screen_seaice_MODIS(modis_loaded_CF0,opts);
     
     modis_loaded_CF0_2 =  load(modis_Nd_file_CF0_2,'CTH','W_time3','Cloud_Fraction_Liquid','Droplet_Number_Concentration_37',...
         'Cloud_Fraction_NoOpt','Cloud_Fraction_Combined','Cloud_Fraction_Undetermined',...
         'gcm_Plat2D_AMSRE','gcm_Plon2D_AMSRE','modisyear_timeseries3','daynum_timeseries3'); 
     opts.years_multi=unique(modis_loaded_CF0_2.modisyear_timeseries3);
     [modis_loaded_CF0_2]=screen_seaice_MODIS(modis_loaded_CF0_2,opts);
     
 end



%CALIOSP CF data - load using script :- read_calipso_monthly_night_IPSL.m for average values
    %first. This puts the data into these fields :- 
    % cllcalipso_monthly_AVERAGE, clmcalipso_monthly_AVERAGE, clhcalipso_monthly_AVERAGE
    % Is average of day and night values for now (have separate ones too).
if ical_data==1
    read_calipso_monthly_night_IPSL
    [out, time_out_CAL, time_inds_CAL, dtime_match] = get_time_range_of_array(array_in,gcm_time_matlab_CALIPSO_monthly,time_choice,dim);
end

if iceres_data==1
    ceres_file = '/home/disk/eos15/d.grosvenor/eos8/CERES/ACSIS/CERES_EBAF-TOA_Ed4.0_Subset_200903-201004.nc'
    nc_ceres = netcdf(ceres_file);
    lon_ceres = nc_ceres{'lon'}(:);
    lat_ceres = nc_ceres{'lat'}(:);    
    time_ceres = nc_ceres{'time'}(:); %days since 2000-03-01
    SW_TOA_ceres = nc_ceres{'toa_sw_all_mon'}(:);
    %solar_mon:long_name = "Incoming Solar Flux, Monthly Means" ;
    SW_TOA_in_ceres = nc_ceres{'solar_mon'}(:);
    cf_ceres = nc_ceres{'cldarea_total_daynight_mon'}(:);
    tau_ceres = nc_ceres{'cldtau_total_day_mon'}(:);
    
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


iread_MAC=0;
if iread_MAC==1
    
   [MAC_lwp_out_2009,MAC_tlwp_out_2009,MAC_time_2009,MAC_lat,MAC_lon,MAC_lat2d,MAC_lon2d] = read_MAC_daily('/home/disk/eos5/d.grosvenor/MAC_LWP/y2009');
   [MAC_lwp_out_2010,MAC_tlwp_out_2010,MAC_time_2010,MAC_lat,MAC_lon,MAC_lat2d,MAC_lon2d] = read_MAC_daily('/home/disk/eos5/d.grosvenor/MAC_LWP/y2010');
   
   MAC_lwp_out_2009_2010 = cat(1,MAC_lwp_out_2009,MAC_lwp_out_2010);
   MAC_tlwp_out_2009_2010 = cat(1,MAC_tlwp_out_2009,MAC_tlwp_out_2010);
   MAC_time_2009_2010 = cat(1,MAC_time_2009,MAC_time_2010);
    
   [out, time_out_MAC, time_inds_MAC, dtime_match] = get_time_range_of_array(array_in,MAC_time_2009_2010,time_choice,dim);
   
   MAC_LWP_time_period_mean = meanNoNan( MAC_lwp_out_2009_2010(time_inds_MAC,:,:),1);
   MAC_TLWP_time_period_mean = meanNoNan( MAC_tlwp_out_2009_2010(time_inds_MAC,:,:),1);           
   
end



%% SW direct and indirect forcing 
isave_plot_global_forcing_draft=1; %for paper draft

isave_plot_global_diff=0;
isave_plot_global_forcing=0;

% Will just consider clouds in the PI run for now
% consider coarse graining to help make sure met is the same?

%Clean means no direct effects (background clean aerosol when out of
%clouds)
% Direct effect is the total (indirect + direct) minus SW_clean (just
% indirect effect, but no direct), i.e. direct = indirect+direct - indirect

 var_UM = 'SW_down_clean_surf'; 
 
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
var_UM = 'SW_down_clean_clear_surf'; 
    
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
     
        
%    dat_modis = meanNoNan(low_CF_PD_ALL,3) - meanNoNan(low_CF_PI_ALL,3);
    
    var_UM = 'Direct forcing (W m^{-2})';
    dat_modis = meanNoNan(direct_ALL,3);    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-10 10]);
    
    % Save
    %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/model_vs_MODIS_Nd_plots/' modis_day_of_year_str '_global_' titlenam_driver];
%    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];

    if isave_plot_global_forcing_draft==1; %for paper draft==1
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %close(gcf);
    end
    
%% Plot indirect
isave_plot_global_forcing_draft=1; %for paper draft

isave_plot_global_diff=0;

    var_UM = 'Indirect forcing (W m^{-2})';
    dat_modis = meanNoNan(indirect_ALL,3);  
    
        %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%    increase_font_size_map_figures
     caxis([-10 10]);
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    
    %UM_ACSIS_run_plot_box_commands  %Plot the box for the regional stats, etc.
    
    % Save  
    if isave_plot_global_forcing_draft==1
%        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end
    
%%     
    
    isave_Dan_Jones=0
    if isave_Dan_Jones==1

        %Save the forcing for Dan Jones
        time_range = 24*(datenum(time_choice.time_range - datenum('01-Jan-1970'))); %Hrs since 1st Jan 1970
        indirect_forcing_time_mean = dat_modis;
        lat2D = gcm_Plat2D_UM;
        lon2D = gcm_Plon2D_UM;

        save_file_Dan_Jones = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/indirect_forcing_timemean.mat'];
        save(save_file_Dan_Jones,'indirect_forcing_time_mean','lat2D','lon2D','time_range','-V7.3');
        mat2nc_Dan(save_file_Dan_Jones,[save_file_Dan_Jones '.nc']);

    end
    
    isave_KH=0;
    if isave_KH==1
        %save forcing for cyclone compositing (Kevin Hodges)
        
        %Round to the nearest hour for consistency with the MSLP timestamps
        [Y,M,D,HH,MM,SS] = datevec(time_out);
        time_out_round = datenum(Y,M,D,HH,0,0);
        time = 24 * (time_out_round - datenum('01-Jan-1970'));        

        
        lat = gcm_Plat2D_UM(:,1);
        long = gcm_Plon2D_UM(1,:);
        
        %Change order from lat,lon,time to time,lat,lon
        indirect_surface_forcing = permute(indirect_ALL,[3 1 2]);
        
        save_file_KH = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/indirect_forcing_ALL_KH.mat'];
        save(save_file_KH,'indirect_surface_forcing','time','lat','long','-V7.3');
        mat2nc_Dan(save_file_KH,[save_file_KH '.nc']);
        eval(['!ncatted -a units,''time'',c,c,''hours since 1970-1-1 00:00:00'' ' save_file_KH '.nc']);
        
        %Rename the dim_180 etc dimension names to proper names
        nt_str = num2str(length(time));
        nlat_str = num2str(length(lat));
        nlon_str = num2str(length(long));
        
        eval_str=['!ncrename -d dim_' nt_str ',time ' save_file_KH '.nc']; eval(eval_str);
        eval_str=['!ncrename -d dim_' nlat_str ',lat ' save_file_KH '.nc']; eval(eval_str);
        eval_str=['!ncrename -d dim_' nlon_str ',long ' save_file_KH '.nc']; eval(eval_str);        
        
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        close(gcf);
    end


%% Make a timeseries
if close_all_figs==1
    close all
    clear gca
end
isave_plot_global_forcing_timser=0;

    %(restricting to the region of interest)
    [indirect_ALL_timser] = UM_make_regional_timeseries(indirect_ALL,nT_main,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM);
    
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        close(gcf);
    end

%% Change in SW out TOA and SW out TOA vs CERES
isave_plot_CERES=1; %For draft paper

if close_all_figs==1
    close all
    clear gca
end

 switch load_type
        case 'cam6'
            
     otherwise

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
        
    SW_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    SW_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    
    SW_TOA_model_save = um_data;
end

 end


%% Plot cloud fractions - low cloud
isave_plot_global_lowCF=0;
isave_plot_CF = 1; %Plot the 3 panel comparisons between UM and CALIPSO as used in paper draft.

if close_all_figs==1
    close all
    clear gca    
end



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
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-0.1 0.1]);
    
    dat_modis = meanNoNan(low_CF_PD_ALL,3); var_UM = 'Low CF in PD)';
    %    dat_modis = meanNoNan(low_CF_PI_ALL,3);
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    %lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    lb_map = lbmap(32,'brownblue'); 
    lb_map = lb_map(1:16,:);
    colormap(flipdim(lb_map,1));
    %colormap(lb_map);
    caxis([0 1]);
    
    
    
    % Save
    if isave_plot_global_lowCF==1
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        close(gcf);
    end

    switch load_type
        case 'cam6'
            
        otherwise
            
            cloud_input2='CALIPSO';
            [low_calipsoCF_PD_ALL,nT_CAL] = UM_ACSIS_global_low_cloud_FUNC('low',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
            
            cloud_input2='MODIS';
            [low_modisCF_PD_ALL,nT_modis] = UM_ACSIS_global_low_cloud_FUNC('low',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
            
            
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
    
    lowCF_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    lowCF_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    
end
    
    
%% mid cloud
isave_plot_global_midCF=0;
isave_plot_CF = 1;

  switch load_type
        case 'cam6'
            
        otherwise

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
    

% Plot cloud fractions - mid cloud

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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        close(gcf);
    end
    
    
    cloud_input2='CALIPSO';
    [mid_calipsoCF_PD_ALL,nT_CAL] = UM_ACSIS_global_low_cloud_FUNC('mid',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
    
    cloud_input2='MODIS';
    [mid_modisCF_PD_ALL,nT_modis] = UM_ACSIS_global_low_cloud_FUNC('mid',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
    
    
    
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
        midCF_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
        midCF_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.

    end
    
    
  end
    
%% High cloud        
isave_plot_global_highCF=0;
isave_plot_CF = 1;

  switch load_type
        case 'cam6'
            
        otherwise

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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        close(gcf);
    end
    
    
    cloud_input2='CALIPSO';
    [high_calipsoCF_PD_ALL,nT_CAL] = UM_ACSIS_global_low_cloud_FUNC('high',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
    
    cloud_input2='MODIS';
    [high_modisCF_PD_ALL,nT_modis] = UM_ACSIS_global_low_cloud_FUNC('high',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
        
    
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
        highCF_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
        highCF_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.

    end
    
  end

    
%% Total cloud   
isave_plot_global_totalCF=0;
icalc_total=1;

  switch load_type
        case 'cam6'
            
        otherwise

    switch cloud_input
        case 'UM'
            var_UM = 'total_max_random_cloud_amount'; %Can get this from stash 9-217 (is added in Python script now).
            icosp_mask=0;                                    
            icalc_total=1;
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
        
        um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        [total_CF_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        total_CF_PI_ALL(total_CF_PI_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud??        



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
            saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
            close(gcf);
        end
        
        
    cloud_input2='CALIPSO';
    [total_calipsoCF_PD_ALL,nT_CAL] = UM_ACSIS_global_low_cloud_FUNC('total',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
    
    cloud_input2='MODIS';
    [total_modisCF_PD_ALL,nT_modis] = UM_ACSIS_global_low_cloud_FUNC('total',cloud_input2,um_case_PD,load_type,time_inds,gcm_Plat2D_UM);
            

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
    
  end

  
%% Plot MODIS COSP optical depth etc

 switch load_type
        case 'cam6'
            
     otherwise
         
         var_UM_mask = ['modis_misr_issp_cloud_amount_mask'];
         um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
         dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
         dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
         [modis_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
         
         var_UM_mask = ['modis_misr_issp_cloud_amount_mask'];
         um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
         dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
         dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
         [modis_mask_PI,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
         
         %Some the quantities need to be weighted by the MODIS total CF.
         var_UM = 'modis_total_cloud_amount'
         um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
         dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
         dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
         [totCF_COSP_MODIS_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %
         totCF_COSP_MODIS_PD_ALL(modis_mask==0) = NaN;
         
         var_UM = 'modis_liquid_cloud_fraction'; %Some of the quantities are multiplied by this (e.g. liquid tau), so need to divide by this
         % to get the in-cloud values.
         um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
         dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
         dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
         [liqCF_COSP_MODIS_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %
         liqCF_COSP_MODIS_PD_ALL(modis_mask==0) = NaN;
         
         iload_PI_MODIS=0;
         if iload_PI_MODIS==1
             var_UM = 'modis_liquid_cloud_fraction'; %Some of the quantities are multiplied by this (e.g. liquid tau), so need to divide by this
             % to get the in-cloud values.
             um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
             dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
             dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
             [liqCF_COSP_MODIS_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %
             liqCF_COSP_MODIS_PI_ALL(modis_mask_PI==0) = NaN;
             
         end
   

    clear gca
    
    var_UM = 'modis_liq_tau'; %N.B. - this is weighted (multplied) by the liquid cloud fraction (2-452; modis_liquid_cloud_fraction)
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [modis_cosp_tau_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    modis_cosp_tau_PD = modis_cosp_tau_PD ./ liqCF_COSP_MODIS_PD_ALL; %Divide by liquid CF to get in-cloud tau
    
    var_UM = 'modis_liq_reff'; %N.B. - this is weighted (multplied) by the liquid cloud fraction (2-452; modis_liquid_cloud_fraction)
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [modis_cosp_reff_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    modis_cosp_reff_PD = modis_cosp_reff_PD ./ liqCF_COSP_MODIS_PD_ALL; %Divide by liquid CF to get in-cloud reff
      
    
    var_UM = 'modis_CTP'; %N.B. - this is weighted (multplied) by the total cloud fraction (2-451; modis_total_cloud_amount)
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [modis_cosp_CTP_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    modis_cosp_CTP_PD = modis_cosp_CTP_PD ./ totCF_COSP_MODIS_PD_ALL; %Divide by COSP MODIS total CF to get in-cloud CTP

    
    var_UM = 'Max_Cloud_Height'; %Max height of cloud with LWC>=0.05 g/kg (in metres); my diag, not COSP
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [max_cloud_height_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    max_cloud_height_PD(max_cloud_height_PD<-0.99)=NaN; %When there is no cloud the value was set to -1
    
    var_UM = 'Max_Cloud_Height_in_cloud_LWC'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [max_cloud_height_in_cloud_LWC_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    max_cloud_height_in_cloud_LWC_PD(max_cloud_height_in_cloud_LWC_PD<-0.99)=NaN; %When there is no cloud the value was set to -1    
    
    var_UM = 'Max_Cloud_Height_in_cloud_LWC_IWC'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [max_cloud_height_in_cloud_LWC_IWC_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    max_cloud_height_in_cloud_LWC_IWC_PD(max_cloud_height_in_cloud_LWC_IWC_PD<-0.99)=NaN; %When there is no cloud the value was set to -1 
    
    iload_MODIS_PI=0;
    if iload_MODIS_PI==1
        var_UM = 'Max_Cloud_Height_in_cloud_LWC_IWC'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
        um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        [max_cloud_height_in_cloud_LWC_IWC_PI,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        max_cloud_height_in_cloud_LWC_IWC_PI(max_cloud_height_in_cloud_LWC_IWC_PI<-0.99)=NaN; %When there is no cloud the value was set to -1
    end
    
    
    switch um_case_PD
        case 'u-bi194 -ignore now'
            var_UM = 'Max_Cloud_Height_in_cloud_LWC_CFice'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
            um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
            dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
            dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
            [max_cloud_height_in_cloud_LWC_CFice_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
            max_cloud_height_in_cloud_LWC_CFice_PD(max_cloud_height_in_cloud_LWC_CFice_PD<-0.99)=NaN; %When there is no cloud the value was set to -1
            
            var_UM = 'Max_CFice'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
            um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
            dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
            dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
            [Max_CFice_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
            Max_CFice_PD(Max_CFice_PD<-0.99)=NaN; %When there is no cloud the value was set to -1
            
            
            dat_modis = max_cloud_height_in_cloud_LWC_CFice_PD(:,:,1); var_UM = 't=1  max_cloud_height_in_cloud_LWC_CFice_PD (Present Day)';            
            %run plotting script
            figure
            UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
            lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
            %caxis([0 1]);
            
            dat_modis = Max_CFice_PD(:,:,1); var_UM = 't=1  Max_CFice_PD (Present Day)';
            %run plotting script
            figure
            UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
            lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
            %caxis([0 1]);            
    end
       
     
    
    %dat_modis = meanNoNan(modis_cosp_tau_PD,3); var_UM = 'Mean MODIS COSP liquid optical depth (Present Day)';
    dat_modis = modis_cosp_tau_PD(:,:,1); var_UM = 't=1 MODIS COSP liquid optical depth (Present Day)';
%    dat_modis = meanNoNan(SW_TOA_PD_ALL,3);
%    dat_modis = meanNoNan(SW_TOA_PI_ALL,3);

    %dat_modis = low_modisCF_PD_ALL(:,:,1); var_UM = 't=1 MODIS COSP low CF (Present Day)';
    dat_modis  = meanNoNan(low_modisCF_PD_ALL(:,:,:),3); var_UM = 'Time mean MODIS COSP low CF (Present Day)';
    
    %dat_modis = mid_modisCF_PD_ALL(:,:,1); var_UM = 't=1 MODIS COSP mid CF (Present Day)';
    
    %dat_modis = high_modisCF_PD_ALL(:,:,1); var_UM = 't=1 MODIS COSP high CF (Present Day)';
    
    low_mid_MODIS = low_modisCF_PD_ALL+mid_modisCF_PD_ALL;
    %dat_modis  = low_mid_MODIS(:,:,1); var_UM = 't=1 MODIS COSP low+mid CF (Present Day)';
    %dat_modis  = meanNoNan(low_mid_MODIS(:,:,:),3); var_UM = 'Time mean MODIS COSP low+mid CF (Present Day)';
    
    low_mid_CAL = low_calipsoCF_PD_ALL+mid_calipsoCF_PD_ALL;
    %dat_modis = meanNoNan(low_mid_CAL,3); var_UM = 'Time mean CALIPSO COSP low+mid CF (Present Day)';    
    %dat_modis = meanNoNan(low_calipsoCF_PD_ALL,3); var_UM = 'Time mean CALIPSO COSP low CF (Present Day)';    
    flow_CAL = low_calipsoCF_PD_ALL ./ low_mid_CAL;
    flow_CAL(low_mid_CAL<0.01)=0;
    %dat_modis = flow_CAL; var_UM = 't=1 CALIPSO COSP fraction low CF to low+mid (Present Day)';
    
    %dat_modis = flow_CAL(:,:,1) .* low_mid_MODIS(:,:,1); var_UM = 't=1 MODIS COSP scaled low CF (Present Day)';
    %dat_modis = meanNoNan(flow_CAL .* low_mid_MODIS , 3); var_UM = 'Time mean MODIS COSP scaled low CF (Present Day)';

    
    CTP_temp = modis_cosp_CTP_PD;
    %CTP_temp(CTP_temp<680e2)=NaN;    
    %dat_modis = CTP_temp(:,:,1)/1e2; var_UM = 't=1 MODIS COSP CTP (Present Day; hPa)';
        
    %dat_modis = low_calipsoCF_PD_ALL(:,:,1); var_UM = 't=1 CALIPSO COSP low CF (Present Day)';
    %dat_modis = mid_calipsoCF_PD_ALL(:,:,1); var_UM = 't=1 CALIPSO COSP mid CF (Present Day)';
    
    %dat_modis = liqCF_COSP_MODIS_PD_ALL(:,:,1); var_UM = 't=1 COSP MODIS liq CF (Present Day)';  
    %dat_modis = meanNoNan(liqCF_COSP_MODIS_PD_ALL(:,:,:) , 3); var_UM = 'Time mean COSP MODIS liq CF (Present Day)';
    

    %dat_modis = totCF_COSP_MODIS_PD_ALL(:,:,1); var_UM = 't=1 COSP MODIS total CF (Present Day)';  
    %dat_modis = meanNoNan(totCF_COSP_MODIS_PD_ALL(:,:,:) , 3); var_UM = 'Time mean COSP MODIS total CF (Present Day)';
        

    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 1]);   

 end
 
 %% Max model cloud height
  dat_modis = max_cloud_height_PD(:,:,1)/1e3; var_UM = 't=1 Max model cloud height (0.05 g/kg threshold; Present Day)';
  dat_modis = meanNoNan(max_cloud_height_PD(:,:,:),3)/1e3; var_UM = 'Time average max model cloud height (0.05 g/kg threshold; Present Day)';    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 4]);  
    
 %% Max model cloud height using in-cloud LWC and IWC
  dat_modis = max_cloud_height_in_cloud_LWC_IWC_PD(:,:,1)/1e3; var_UM = 't=1 Max model cloud height (0.05 g/kg IN-CLOUD threshold; Present Day)';
  dat_modis = meanNoNan(max_cloud_height_in_cloud_LWC_IWC_PD(:,:,:),3)/1e3; var_UM = 'Time average max model cloud height (0.05 g/kg IN-CLOUD threshold; Present Day)';    
  %dat_modis = meanNoNan(max_cloud_height_in_cloud_LWC_CFice_PD(:,:,:),3)/1e3; var_UM = 'Time average max model cloud height (0.05 g/kg IN-CLOUD threshold; Present Day)';      
   %CFice method doesn't work very well
   
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    %lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    lb_map = lbmap(32,'brownblue'); 
    lb_map = lb_map(1:16,:);
    colormap(flipdim(lb_map,1));    
    colormap(lb_map);
    caxis([0 4]);      
    
    titlenam_driver = 'CTH_mean';
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
    
 %% CTP for CTH close to 3.2 km
  UM_CTP = modis_cosp_CTP_PD;
  target=3.2e3;
  tol=0.2e3;
  inan = find(max_cloud_height_PD < target-tol | max_cloud_height_PD > target+tol | isnan(max_cloud_height_PD)==1 );  
  UM_CTP(inan)=NaN;
  dat_modis = UM_CTP(:,:,1)/1e2; var_UM = 't=1 Mean MODIS COSP CTP for height match';
  dat_modis = meanNoNan(UM_CTP(:,:,:),3)/1e2; var_UM = 'Time average MODIS COSP CTP for height match';    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([300 1100]);     
  

    iRyanE_predictors=0;
%% Predictors for Ryan E            
    if iRyanE_predictors==1
    
    var_UM = 'LTS'; % Predictor variables for Ryan E.
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [LTS_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); 
    
    var_UM = 'omega700_multi'; % Predictor variables for Ryan E.
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [omega700_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    var_UM = 'qv700_multi'; % in kg/kg
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [qv700_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    var_UM = 'surface_temp'; % Predictor variables for Ryan E.
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [surface_temp_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    var_UM = 'U10'; % Predictor variables for Ryan E.
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global_UV = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [U10_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global_UV,time_inds,load_type,gcm_Plat2D_UM);
            
    var_UM = 'V10'; % Predictor variables for Ryan E.
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [V10_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    
    
    end
    
    
    
%% Droplet concs   
isave_plot_global_Nd = 1;
isave_plot_Nd_bias=0; %Used for the draft paper (3 panel bias plot)
icoarse_grain=0;


if close_all_figs==1
    close all
    clear gca    
end

 Nd_type = 'lwc weighted';
 Nd_type = 'lwc in-cloud weighted'; 
 Nd_type = 'lwc in-cloud weighted ztop'; 
 
    switch Nd_type
        case 'lwc weighted'
            var_UM_Nd = 'Nd_lwc_weighted_UKCA';
        case 'lwc in-cloud weighted'
            var_UM_Nd = 'Nd_lwc_in_cloud_weighted_UKCA';
        case 'lwc in-cloud weighted ztop'
            var_UM_Nd = 'Nd_lwc_in_cloud_weighted_UKCA_ztop';            
        
    end
    
    icosp_mask=0;
    
        
 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM_Nd,run_type,load_type);   
    [Nd_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
%    total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 
    
    if icosp_mask==1
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type);
        [total_CF_mask,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
        total_CF_PD_ALL(total_CF_mask==0) = NaN;
    end
    

    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM_Nd,run_type,load_type);   
    [Nd_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
%    total_CF_PD_ALL(total_CF_PD_ALL>1.001)=NaN; %got some strange spikes in the data for mid-cloud?? 


 
    
    

    
%     dat_modis = meanNoNan(clhcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
%     figure
%     %run plotting script
%     CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%     lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%     caxis([0 1]);
    
    
    %Run subplotting script for UM, MODIS and bias
    UM_ACSIS_global_SUBPLOT_commands_Nd    %saves here too
    
%Change in Nd from PI to PD

    Nd_PI_map = meanNoNan(Nd_PI_ALL,3);
    Nd_PD_map = meanNoNan(Nd_PD_ALL,3);
      
    
    dat_modis = (Nd_PD_map - Nd_PI_map)/1e6; var_UM = 'Change in Nd (PD minus PI; cm^{-3})';    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]);
    
if exist('um_case_CONV_SCAV_ON')
    um_case=um_case_CONV_SCAV_ON; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM_Nd,run_type,load_type);   
    [Nd_CONV_SCAV_ON_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);  %per m3
    Nd_CONV_SCAV_ON_ALL_map = meanNoNan(Nd_CONV_SCAV_ON_ALL,3);
 
    
    dat_modis = (Nd_CONV_SCAV_ON_ALL_map - Nd_PD_map)/1e6; var_UM = 'Change in Nd (Conv scav ON minus all scav OFF; cm^{-3})';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]); 
    
    dat_modis = (Nd_PI_map - Nd_CONV_SCAV_ON_ALL_map )/1e6; var_UM = 'Change in Nd (ALL scav ON minus just conv scav ON; cm^{-3})';
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]); 
    
    end   
    
    % Save

    if isave_plot_global_Nd==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%        close(gcf);
    end
    
    dat_modis = 100 * (Nd_PD_map - Nd_PI_map) ./ Nd_PI_map ; var_UM = 'Percentage change in Nd (PD minus PI)';    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-150 150]);  

    Nd_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    Nd_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    
    
    
    if isave_plot_global_Nd==1
        savename=[savedir_date titlenam_driver];       
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%        close(gcf);
    end
    
    
    %PD Nd map
    dat_modis = (Nd_PD_map)/1e6; var_UM = ['Nd in PD for ' um_case_PD '; (cm^{-3})'];    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(32,'brownblue'); 
    lb_map = lb_map(1:16,:);
    colormap(flipdim(lb_map,1));
    caxis([0 300]);
    
    titlenam_driver = var_UM;
    savename=[savedir_date titlenam_driver];       
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        
        
    
    %PI Nd map
    dat_modis = (Nd_PI_map)/1e6; var_UM = ['Nd in PI for ' um_case_PI '; (cm^{-3})'];    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-100 100]);
    
%% Nd longitudinal transects etc.   
iplot_long_transects=0;
if iplot_long_transects==1

    clear leg_str
    lat_find = -20;
    lon_range = [-120 -70];
    [minval, ilat] = min(abs(gcm_Plat2D_UM(:,1) - lat_find));
    ilons = find(gcm_Plon2D_UM(1,:) > lon_range(1) & gcm_Plon2D_UM(1,:) <= lon_range(2));
    lons = gcm_Plon2D_UM(1,ilons);
    Nd_trans_PI = Nd_PI_map(ilat,ilons)/1e6; 
    Nd_trans_PD = Nd_PD_map(ilat,ilons)/1e6; 

    figure
    clear leg_str
    leg_str{1} = 'Scavenging OFF';
    leg_str{2} = 'Scavenging ON';
    plot(lons,Nd_trans_PI,'r','linewidth',3); hold on; grid on
    plot(lons,Nd_trans_PD,'r--','linewidth',3);    
    set(gca,'fontsize',18);
    legend(leg_str,'location','northwest');
    ylabel('Droplet Concentration (cm^{-3})');
    xlabel('Longitude (degrees)');
    title('Transect at lat=20^{o}S');
    %increase_font_size_map_figures        
    
    titlenam_driver='Nd_transect_plot';
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
    
    cth_mean_PD = meanNoNan(max_cloud_height_in_cloud_LWC_IWC_PD,3);
    cth_mean_PI = meanNoNan(max_cloud_height_in_cloud_LWC_IWC_PI,3);
    
    cth_trans_PI = cth_mean_PI(ilat,ilons); 
    cth_trans_PD = cth_mean_PD(ilat,ilons); 

    figure
    clear leg_str
    leg_str{1} = 'Scavenging OFF';
    leg_str{2} = 'Scavenging ON';
    plot(lons,cth_trans_PI,'r','linewidth',3); hold on; grid on
    plot(lons,cth_trans_PD,'r--','linewidth',3);    
    set(gca,'fontsize',18);
    legend(leg_str);
    ylabel('Cloud Top Height (m)');
    xlabel('Longitude (degrees)');
    title('Transect at lat=20^{o}S');
    
    titlenam_driver='CTH_transect_plot';
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
         
    
    figure
    clear leg_str
    leg_str{1} = 'Scavenging OFF';
    leg_str{2} = 'Scavenging ON';
    plot(lons,1e6*Nd_trans_PI .* cth_trans_PI,'r','linewidth',3); hold on; grid on
    plot(lons,1e6*Nd_trans_PD .* cth_trans_PD,'r--','linewidth',3);    
    set(gca,'fontsize',18);
    legend(leg_str);
    ylabel('Nd * CTH (m^{-2})');
    xlabel('Longitude (degrees)');
    title('Transect at lat=20^{o}S');  
    
    figure
    clear leg_str
    leg_str{1} = 'Scavenging OFF';
    leg_str{2} = 'Scavenging ON';
    leg_str{3} = 'Scavenging OFF, Nd_{CTH}'; 
    leg_str{4} = 'Scavenging ON, Nd_{CTH}';
    plot(lons,Nd_trans_PI,'r','linewidth',3); hold on; grid on
    plot(lons,Nd_trans_PD,'r--','linewidth',3);    
    plot(lons,2.25e5 ./ cth_trans_PI,'b','linewidth',3); hold on; grid on
    plot(lons,2.00e5 ./ cth_trans_PD,'b--','linewidth',3);    
    set(gca,'fontsize',18);
    
    legend(leg_str,'location','northwest');
    ylabel('Nd_{CTH} (cm^{-3})');
    xlabel('Longitude (degrees)');
    title('Transect at lat=20^{o}S');  
    
   
    

 
    
    figure
     clear leg_str   
     leg_str{1} = 'Scavenging OFF';
     leg_str{2} = 'Scavenging OFF, Nd_{CTH}';
    plot(lons,Nd_trans_PI,'r','linewidth',3); hold on; grid on
    %plot(lons,Nd_trans_PD,'r--','linewidth',3);    
    plot(lons,2.25e5 ./ cth_trans_PI,'b','linewidth',3); hold on; grid on
    %plot(lons,2.25e5 ./ cth_trans_PD,'b--','linewidth',3);    
    set(gca,'fontsize',18);        
    legend(leg_str,'location','northwest');
    ylabel('Nd_{CTH} (cm^{-3})');
    xlabel('Longitude (degrees)');
    title('Transect at lat=20^{o}S');  
    
     titlenam_driver='Nd_estimate_scav_transect_plot';
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
    

figure
    clear leg_str 
    leg_str{1} = 'Scavenging ON';
    leg_str{2} = 'Scavenging ON, Nd_{CTH}';
    %plot(lons,Nd_trans_PI,'r','linewidth',3); hold on; grid on
    plot(lons,Nd_trans_PD,'r--','linewidth',3);
    hold on; grid on
    %plot(lons,2e5 ./ cth_trans_PI,'b','linewidth',3); hold on; grid on
    plot(lons,2.25e5 ./ cth_trans_PD,'b--','linewidth',3);    
    set(gca,'fontsize',18);    
    legend(leg_str,'location','northwest');
    ylabel('Nd_{CTH} (cm^{-3})');
    xlabel('Longitude (degrees)');
    title('Transect at lat=20^{o}S');   
    
    titlenam_driver='Nd_estimate_no_scav_transect_plot';
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);  
    
end
    

%% Minghuai Nd trends

iplot_Ming_trends=0;
if iplot_Ming_trends==1
    
var_UM = 'Nd_cf_weighted_UKESM_ztop';
um_case='UKESM/r1i1p1f2_u-bc179/'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

years=[2000:2014];
clear Nd_annual
for iy=1:15
    tind_01 = (iy-1)*12+1;
    tind_02 = tind_01+11;
    for im=1:12
        Nd_annual(iy,:,:) = meanNoNan(dat_global.dat(tind_01:tind_02,:,:),1);
    end
end



 %y2015 minus y2000
    istart=1; iend=15;
    dat_modis = squeeze (Nd_annual(iend,:,:) - Nd_annual(istart,:,:) ) / (iend-istart); 
    var_UM = ['Nd trend between y' num2str(years(istart)) ' and y' num2str(years(iend)) '; (cm^{-3} yr^{-1})'];    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);
        
    istart=2; iend=15;
    dat_modis = squeeze (Nd_annual(iend,:,:) - Nd_annual(istart,:,:) ) / (iend-istart);     
    var_UM = ['Nd trend between y' num2str(years(istart)) ' and y' num2str(years(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);
            
    istart=1; iend=14;
    dat_modis = squeeze (Nd_annual(iend,:,:) - Nd_annual(istart,:,:) ) / (iend-istart);     
    var_UM = ['Nd trend between y' num2str(years(istart)) ' and y' num2str(years(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);    

end    
        
%% Satellite LWP
    
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
    
    
    
%% Plot LWP
isave_plot_global_LWP=0;
isave_plot_LWP_bias=1; %draft paper 3 panel plot

clear gca

% var_UM = 'LWP'; %calculated from Python script using QCL after timestep(0-254)
 var_UM = 'LWP_sec30'; %from section 30-405

 
    um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);  
    [LWP_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);

    
 var_UM = 'LWP_sec30';  %from section 30-405
 
    um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);   
    [LWP_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);   
    
    switch var_UM
        case {'LWP_sec30','LWP_correct_dz'}
           LWP_PD_ALL =  LWP_PD_ALL*1e3; %convert to g/m2
           LWP_PI_ALL =  LWP_PI_ALL*1e3; %convert to g/m2           
    end
    
    dat_modis = meanNoNan(LWP_PI_ALL,3);
    
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%    caxis([0 1]);
    
    % Save
    if isave_plot_global_LWP==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        close(gcf);
    end
    
%% RWP part.    
    
    
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
    
    icoarse_grain=0;

    iload_conv_lwp_rwp=1;
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


        icoarse_grain=0;
    
%% Convective LWP fraction from model.    
    %Plot the fraction of LWP that comes from the convection scheme
    TLWP_LS = meanNoNan(LWP_PD_ALL + RWP_PD_ALL,3);
    TLWP_CONV = meanNoNan(Conv_LWP_PD_ALL + Conv_RWP_PD_ALL,3);
    dat_modis = TLWP_CONV ./ (TLWP_LS + TLWP_CONV);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    icoarse_grain=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    icoarse_grain=1;
    caxis([0 1]);
    title('Fraction of LWP+RWP from convection scheme');
    
    isave_plot_conv_LWP_fraction=1;
     if isave_plot_conv_LWP_fraction==1
%        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_LWP_bias_3panel_vs_satellite'];

        titlenam_driver=['conv_TLWP_fraction'];
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        
%         iappend=0;
%         for i=1:length(region_biases)
%             latex_newcommand_from_structure(region_biases{i},[region_biases{i}.region_shortname var_Sc],savename,iappend);
%             iappend=1;
%         end
        
     end
    
    %Plot the fraction of LWP that comes from the convection scheme
    LWP_LS = meanNoNan(LWP_PD_ALL,3);
    LWP_CONV = meanNoNan(Conv_LWP_PD_ALL,3);
    dat_modis = LWP_CONV ./ (LWP_LS + LWP_CONV);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    icoarse_grain=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    icoarse_grain=1;
    caxis([0 1]);
    title('Fraction of LWP from convection scheme');
    
        %Plot the fraction of LWP that comes from the convection scheme    
    dat_modis = LWP_CONV ./ (TLWP_LS + TLWP_CONV);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    icoarse_grain=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    icoarse_grain=1;
    caxis([0 1]);
    title('Fraction of LWP from convection scheme to TLWP');
    
    %Plot the fraction of LWP that comes from the convection scheme    
    RWP_CONV = meanNoNan(Conv_RWP_PD_ALL,3);
    dat_modis = RWP_CONV ./ (TLWP_LS + TLWP_CONV);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    icoarse_grain=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    icoarse_grain=1;
    caxis([0 1]);
    title('Fraction of TLWP that is RWP from convection scheme');
    icoarse_grain=0; %reset back

%%     
    %Same for PI
    TLWP_LS_PI = meanNoNan(LWP_PI_ALL + RWP_PI_ALL,3);
    TLWP_CONV_PI = meanNoNan(Conv_LWP_PI_ALL + Conv_RWP_PI_ALL,3);
    dat_modis = TLWP_CONV_PI ./ (TLWP_LS_PI + TLWP_CONV_PI);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 1]);
    title('Fraction of LWP+RWP from convection scheme for PI');
    
%% LWP to TLWP fraction from model inc convection scheme.    
    %Plot the fraction of LWP that comes from the convection scheme   
    LWP_LS_CONV = LWP_PD_ALL + Conv_LWP_PD_ALL;
    TLWP_LS_CONV = LWP_LS_CONV + RWP_PD_ALL + Conv_RWP_PD_ALL;
    inan=find(TLWP_LS_CONV<0.1);
    TLWP_LS_CONV2 = TLWP_LS_CONV;
    TLWP_LS_CONV2(inan) = NaN;
        
    dat_modis =  meanNoNan(LWP_LS_CONV,3) ./ meanNoNan(TLWP_LS_CONV2,3);    
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    icoarse_grain=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    icoarse_grain=1;
    caxis([0 1]);
    title('Fraction of LWP to LWP+RWP inc. convection scheme');

    %% LWP to TLWP fraction from model WITHOUT convection scheme.    
    %Plot the fraction of LWP that comes from the convection scheme   
    %LWP_LS_ALL = LWP_PD_ALL;
    TLWP_LS_no_CONV = LWP_PD_ALL + RWP_PD_ALL;
    inan=find(TLWP_LS_no_CONV<0.1);
    TLWP_LS_no_CONV2 = TLWP_LS_no_CONV;
    TLWP_LS_no_CONV2(inan) = NaN;
        
    dat_modis =  meanNoNan(LWP_PD_ALL,3) ./ meanNoNan(TLWP_LS_no_CONV2,3);    
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    icoarse_grain=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global    
    caxis([0 1]);
    title('Fraction of LWP to LWP+RWP WITHOUT convection scheme');

    
   

end




%% LWP to TWLP ratio bias plot and LWP comparison
isave_plot_LWP_bias_Sc=1;

    %Calculate and plot the ratio first (need this for the contour too).
    ioverride_LWP_Sc_choices=1;
    plot_var_Sc = 'Ratio of LWP to TLWP';
    lwp_day_night = 'average'; %Fairly weak dirunal cycle.
    iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
    irestrict_using_LWP_fraction=0; %whether to restrict to boxes with a high LWP to TLWP ratio
    min_LWP_ratio=0.99;
    iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
    ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
    %mid or high cloudedit plo
    % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1
    cloud_scene_selection = 'none'; model_CTH_thresh = 3.2e3; %CTH threshold in km    
    irestrict_model_CF = 0;  %whether to restric the model CFs to > 80%
    imodis_cf = 0;
    modis_CTH_screening = 'none';
    LWP_sat = 'AMSR-E';
    fsize_latlon = 14;    
    ibias_contour=0;
    bias_type = 'Percentage'; %overrides to absolute if plotting ratio of LWP to TLWP
    min_LWP_ratio=0.99;
    
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
    LWP_ratio_AMSRE = sat_data; %save the LWP to TLWP ratio from the satellite for the contours in the next plot.
    LWP_ratio_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    LWP_ratio_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.

    
    
    
    
    
    %Plot the PD minus PI change of different LWPs and RWPs to check for convective
    %changes  
    
    % LS LWP change
    dat_modis = meanNoNan(LWP_PD_ALL,3) - meanNoNan(LWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-10 10]);
    title('\DeltaLWP (PD minus PI, g m^{-2})'); 
    
    % Percent  LS LWP change
    dat_modis = 100* ( meanNoNan(LWP_PD_ALL,3) - meanNoNan(LWP_PI_ALL,3) ) ./ meanNoNan(LWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage LWP change (PD vs PI)'); 
    
    % LS RWP change
    dat_modis = meanNoNan(RWP_PD_ALL,3) - meanNoNan(RWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-10 10]);
    title('\DeltaRWP (PD minus PI, g m^{-2})'); 
    
    % Percent  LS RWP change
    dat_modis = 100* ( meanNoNan(RWP_PD_ALL,3) - meanNoNan(RWP_PI_ALL,3) ) ./ meanNoNan(RWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage RWP change (PD vs PI)');
    
    % in-cloud RWP
    cf_min = 0.01;
    RWP_PI_ALL_ic = RWP_PI_ALL ./ low_CF_PI_ALL;
    RWP_PD_ALL_ic = RWP_PD_ALL ./ low_CF_PD_ALL;    
    RWP_PI_ALL_ic(low_CF_PI_ALL<cf_min)=NaN; RWP_PD_ALL_ic(low_CF_PD_ALL<cf_min)=NaN; %regions with no cloud - set W to NaN since have no cf to divide by

    % Plot percent  LS RWP change
    dat_modis = 100* ( meanNoNan(RWP_PD_ALL_ic,3) - meanNoNan(RWP_PI_ALL_ic,3) ) ./ meanNoNan(RWP_PI_ALL_ic,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage RWP_{in-cloud} change (PD vs PI)');
    
    
    
  if iload_conv_lwp_rwp==1
      
      % Conv LWP change
    dat_modis = meanNoNan(Conv_LWP_PD_ALL,3) - meanNoNan(Conv_LWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-10 10]);
    title('\DeltaLWP_{conv} (PD minus PI, g m^{-2})'); 
    
    % Percent  Conv LWP change
    dat_modis = 100* ( meanNoNan(Conv_LWP_PD_ALL,3) - meanNoNan(Conv_LWP_PI_ALL,3) ) ./ meanNoNan(Conv_LWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage LWP_{conv} change (PD vs PI)'); 
    
    % Conv RWP change
    dat_modis = meanNoNan(Conv_RWP_PD_ALL,3) - meanNoNan(Conv_RWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-10 10]);
    title('\DeltaRWP_{conv} (PD minus PI, g m^{-2})'); 
    
    % Percent Conv RWP change
    dat_modis = 100* ( meanNoNan(Conv_RWP_PD_ALL,3) - meanNoNan(Conv_RWP_PI_ALL,3) ) ./ meanNoNan(Conv_RWP_PI_ALL,3);
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage RWP_{conv} change (PD vs PI)');
    
    %Now the changes in combined LWP (LWP_LS + LWP_Conv) and combined RWP, plus
    %the change of the overall combined LWP and RWP
        
%    TLWP_LS_PI = meanNoNan(LWP_PI_ALL + RWP_PI_ALL,3);
%    TLWP_CONV_PI = meanNoNan(Conv_LWP_PI_ALL + Conv_RWP_PI_ALL,3);

    LWP_combined_PI = meanNoNan(LWP_PI_ALL,3) + meanNoNan(Conv_LWP_PI_ALL,3);
    RWP_combined_PI = meanNoNan(RWP_PI_ALL,3) + meanNoNan(Conv_RWP_PI_ALL,3);   
    LWP_combined_PD = meanNoNan(LWP_PD_ALL,3) + meanNoNan(Conv_LWP_PD_ALL,3);
    RWP_combined_PD = meanNoNan(RWP_PD_ALL,3) + meanNoNan(Conv_RWP_PD_ALL,3);   
    
    % Percent Combined (LS+Conv) LWP change
    dat_modis = 100* ( LWP_combined_PD - LWP_combined_PI ) ./ LWP_combined_PI;
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage LWP_{combined} change (PD vs PI)');
    
    % Percent Combined (LS+Conv) RWP change
    dat_modis = 100* ( RWP_combined_PD - RWP_combined_PI ) ./ RWP_combined_PI;
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage RWP_{combined} change (PD vs PI)');
    
    
    TLWP_combined_PI = LWP_combined_PI + RWP_combined_PI;
    TLWP_combined_PD = LWP_combined_PD + RWP_combined_PD;    
    
    % Percent Combined (LS+Conv) TLWP (RWP+LWP) change
    dat_modis = 100* ( TLWP_combined_PD - TLWP_combined_PI ) ./ TLWP_combined_PI;
    %run plotting script - N.B. - if have run to this just stage and then
    %set irestrict_domain_DRIVER=0 the global plot won't work - have to
    %re-run the lat/lon calculations.
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([-50 50]);
    title('Percentage TLWP_{combined} change (PD vs PI)');
    

    
      
      
  end
    
    
    
    
    

   
%% Do the 3-panel plot of model vs AMSRE    
    %Run subplotting script for UM, CALIPSO and bias
    UM_ACSIS_global_SUBPLOT_commands_LWP
    
 
%% Plot ratio of LWP to TLWP from MAC 
if iread_MAC==1
    
    figure
    MAC_ratio = MAC_LWP_time_period_mean./MAC_TLWP_time_period_mean;
    dat_modis = griddata(MAC_lat2d,MAC_lon2d,MAC_ratio,gcm_Plat2D_UM,gcm_Plon2D_UM);
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 1]);
    title('Ratio of LWP to TLWP from MAC');
    
end
    

%% Now plot the TLWP from the model vs obs (no filtering) with contour of LWP to TLWP ratio
    ioverride_LWP_Sc_choices=1;
    plot_var_Sc = 'LWP';
    lwp_day_night = 'average'; %Fairly weak dirunal cycle.
    iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
    irestrict_using_LWP_fraction=0; %whether to restrict to boxes with a high LWP to TLWP ratio
    %min_LWP_ratio=0.99;
    iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
    ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
    %mid or high cloudedit plo
    % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1
    cloud_scene_selection = 'none';
    irestrict_model_CF = 0;  %whether to restric the model CFs to > 80%
    imodis_cf = 0;
    modis_CTH_screening = 'none';
    LWP_sat = 'AMSR-E';
    fsize_latlon = 14;
    ibias_contour=1;
    bias_type = 'Percentage'; %overrides to absolute if plotting ratio of LWP to TLWP
    
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
    
    LWP_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    LWP_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    
    
%% LWP comparison only for times with high LWP to TLWP fraction
    %Calculate and plot the ratio first (need this for the contour too).
    ioverride_LWP_Sc_choices=1;
    plot_var_Sc = 'LWP';
    lwp_day_night = 'average'; %Fairly weak dirunal cycle.
    iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
    irestrict_using_LWP_fraction=1; %whether to restrict to boxes with a high LWP to TLWP ratio
    min_LWP_ratio=0.99;
    iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
    ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
    %mid or high cloud
    % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1
    cloud_scene_selection = 'none';
    irestrict_model_CF = 0;  %whether to restric the model CFs to > 80%
    imodis_cf = 0;
    modis_CTH_screening = 'none';
    LWP_sat = 'AMSR-E';
    fsize_latlon = 14;    
    ibias_contour=0;
    bias_type = 'Percentage'; %overrides to absolute if plotting ratio of LWP to TLW
    
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
    
    LWP_Sc_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    LWP_Sc_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    
%% In-cloud LWP - Now plot the TLWP from the model vs obs
    ioverride_LWP_Sc_choices=1;
    plot_var_Sc = 'LWPic';
    lwp_day_night = 'average'; %Fairly weak dirunal cycle.
    iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
    irestrict_using_LWP_fraction=0; %whether to restrict to boxes with a high LWP to TLWP ratio
    %min_LWP_ratio=0.99;
    iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
    ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
    %mid or high cloudedit plo
    % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1
    cloud_scene_selection = 'none';
    irestrict_model_CF = 0;  %whether to restric the model CFs to > 80%
    imodis_cf = 0;
    modis_CTH_screening = 'none';
    LWP_sat = 'AMSR-E';
    fsize_latlon = 14;
    ibias_contour=1;
    bias_type = 'Percentage'; %overrides to absolute if plotting ratio of LWP to TLWP
    
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
    
    lwpic_obs_mean = sat_data;
    lwpic_model = um_data; %time average data
    LWPic_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    LWPic_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    %also is thisw saved in the above script
    %lwpic_model_all_times = dat_modis_temp ./ low_calipsoCF_PD_ALL;
    % I.e. the in-cloud LWP for each time individually, rather than using
    % the time mean low CF to divide by
    lwpic_model_all_times_save = lwpic_model_all_times;
    
    
%% In-cloud LWP - LWP comparison only for times with high LWP to TLWP fraction
    %Calculate and plot the ratio first (need this for the contour too).
    ioverride_LWP_Sc_choices=1;
    plot_var_Sc = 'LWPic';
    lwp_day_night = 'average'; %Fairly weak dirunal cycle.
    iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
    irestrict_using_LWP_fraction=1; %whether to restrict to boxes with a high LWP to TLWP ratio
    min_LWP_ratio=0.99;
    iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
    ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
    %mid or high cloud
    % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1
    cloud_scene_selection = 'none';
    irestrict_model_CF = 0;  %whether to restric the model CFs to > 80%
    imodis_cf = 0;
    modis_CTH_screening = 'none';
    LWP_sat = 'AMSR-E';
    fsize_latlon = 14;    
    ibias_contour=0;
    bias_type = 'Percentage'; %overrides to absolute if plotting ratio of LWP to TLW
    
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
    
    lwpic_obs_mean_restrict = sat_data;
    lwpic_model_restrict = um_data; %time average data
    LWPic_Sc_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    LWPic_Sc_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    lwpic_model_all_times_restrict_save = lwpic_model_all_times;

%% Droplet concs filtered for CF>80% etc.
    
    ioverride_LWP_Sc_choices=1;
    plot_var_Sc = 'Nd'; %sets LWP_sat to 'MODIS Nd'
    %lwp_day_night = 'average'; %Fairly weak dirunal cycle.
    %iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
    %irestrict_using_LWP_fraction=0; %whether to restrict to boxes with a high LWP to TLWP ratio
    %min_LWP_ratio=0.99;
    %iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
    %ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
    %mid or high cloudedit plo
    % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1
    cloud_scene_selection = 'Model CTH screening'; model_CTH_thresh = 3.2e3; %CTH threshold in km    
    irestrict_model_CF = 1;  %whether to restrict the model CFs to > 80%
    imodis_cf = 1;
    %modis_CTH_screening = 'none';
    %LWP_sat = 'AMSR-E'; %set to 'MODIS Nd' when plot_var_Sc is set to 'Nd'
    modis_CTH_screening = '3.2km, CF>80';
    fsize_latlon = 14;    
    ibias_contour=0;
    irestrict_using_LWP_fraction=0; %whether to restrict to boxes with a high LWP to TLWP ratio
    
    
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only    
    Nd_prc_bias = prc_bias; %save the bias to do spatial correleations with the SW bias.
    Nd_abs_bias = abs_bias; %save the bias to do spatial correleations with the SW bias.
    Nd_CTH_3pt2km_CF80 = um_data;
    Nd_sat_data_saved = sat_data;
    
    if isave_plot_global_Nd==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%        close(gcf);
    end
    
    
    
    
%% Calculate spatial correlations between the SW TOA biases and the biases in low CF, Nd and LWP
inan=find( isnan(lowCF_prc_bias)==1 | isnan(midCF_prc_bias)==1 | isnan(highCF_prc_bias)==1 | isnan(Nd_prc_bias)==1 | isnan(LWP_prc_bias)==1 | isnan(LWP_Sc_prc_bias)==1 | isnan(SW_prc_bias)==1 );

sat_data = SW_prc_bias; sat_data(inan)=NaN;
iplot_box_on_map=0;

savename=[savedir_date 'SW_bias_lowCF'];
var_Latex = 'SWlowCF_bias';
um_data = lowCF_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'SW_bias_midCF'];
var_Latex = 'SWmidCF_bias';
um_data = midCF_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'SW_bias_highCF'];
var_Latex = 'SWhighCF_bias';
um_data = highCF_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'SW_bias_Nd'];
var_Latex = 'SWNd_bias';
um_data = Nd_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'SW_bias_LWP'];
var_Latex = 'LWP_bias';
um_data = LWP_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'SW_bias_LWP_Sc'];
var_Latex = 'LWP_Sc_bias';
um_data = LWP_Sc_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'SW_bias_LWPic'];
var_Latex = 'LWPic_bias';
um_data = LWPic_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'SW_bias_LWPic_Sc'];
var_Latex = 'LWPic_Sc_bias';
um_data = LWPic_Sc_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

savename=[savedir_date 'CF_LWP_corr'];
var_Latex = 'CF_corr_LWP_bias';
um_data = LWP_prc_bias;
sat_data = lowCF_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region
    
savename=[savedir_date 'CF_LWP_Sc_corr'];
var_Latex = 'CF_corr_LWP_Sc_bias';
um_data = LWP_Sc_prc_bias;
sat_data = lowCF_prc_bias;
DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region

%% Surface rain rates    
if close_all_figs==1
    close all
    clear gca    
end

switch load_type
    case 'cam6'
        
    otherwise
        
        var_UM = 'LS_surf_rain_rate'; %from section 30-405
        
        um_case=um_case_PI; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        [surf_rain_rate_PI_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); %in kg/m2/s
        
        um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        [surf_rain_rate_PD_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM); %in kg/m2/s
        
        % LS surface precip rate in PD
        dat_modis = ( meanNoNan(surf_rain_rate_PD_ALL,3) );
        %run plotting script - N.B. - if have run to this just stage and then
        %set irestrict_domain_DRIVER=0 the global plot won't work - have to
        %re-run the lat/lon calculations.
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        caxis([0 1e-5]);
        title('Surface rain rate (PD, kg m^{-2} hr^{-1})');
        
        % Absolute LS surface precip rate change
        dat_modis = ( meanNoNan(surf_rain_rate_PD_ALL,3) - meanNoNan(surf_rain_rate_PI_ALL,3) );
        %run plotting script - N.B. - if have run to this just stage and then
        %set irestrict_domain_DRIVER=0 the global plot won't work - have to
        %re-run the lat/lon calculations.
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        caxis([-0.02 0.02]);
        title('Change in surface rain rate (PD minus PI, kg m^{-2} hr^{-1})');
        
        var_UM = 'dummy';
        % Percent LS surface precip rate change
        dat_modis = 100* ( meanNoNan(surf_rain_rate_PD_ALL,3) - meanNoNan(surf_rain_rate_PI_ALL,3) ) ./ meanNoNan(surf_rain_rate_PI_ALL,3);
        %run plotting script - N.B. - if have run to this just stage and then
        %set irestrict_domain_DRIVER=0 the global plot won't work - have to
        %re-run the lat/lon calculations.
        figure
        UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
        caxis([-50 50]);
        title('Percentage surface rain rate change (PD vs PI)');
        
end

    
%% SW forcing vs cloud fraction / height, etc.


%% Restrict arrays to the box region(s)
if close_all_figs==1
    close all
    clear gca    
end

 switch load_type
        case 'cam6'
            
     otherwise

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
    
    [iregion_lin,iregion_lin_edges,cf_low]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,low_CF_PI_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_mid]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,mid_CF_PI_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_high]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,high_CF_PI_ALL(:,:,1:end));     
    [iregion_lin,iregion_lin_edges,cf_total]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,total_CF_PI_ALL(:,:,1:end));
        %CF has one extra time dimension for some reason - is offset by a
        %timestep I think - regrid in time?
        
    [iregion_lin,iregion_lin_edges,cf_low_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,low_CF_PD_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_mid_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,mid_CF_PD_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_high_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,high_CF_PD_ALL(:,:,1:end));                 
    [iregion_lin,iregion_lin_edges,cf_total_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,total_CF_PD_ALL(:,:,1:end));
    
    [iregion_lin,iregion_lin_edges,cf_low_av]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,low_CF_av_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_mid_av]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,mid_CF_av_ALL(:,:,1:end)); 
    [iregion_lin,iregion_lin_edges,cf_high_av]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,high_CF_av_ALL(:,:,1:end));             
    
    cf_low_max = max(cf_low,cf_low_PD);
    cf_mid_max = max(cf_mid,cf_mid_PD);
    cf_high_max = max(cf_high,cf_high_PD);    
    
    
    [iregion_lin,iregion_lin_edges,Nd_region]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,Nd_PD_ALL(:,:,1:end));
    [iregion_lin,iregion_lin_edges,liqCF_COSP_region]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,liqCF_COSP_MODIS_PD_ALL(:,:,1:end));
    
 end
    
    
%    forcing_ALL = SW_down_PD_ALL - SW_down_PI_ALL;
    forcing_ALL = indirect_ALL;
    [iregion_lin,iregion_lin_edges,forcing]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,forcing_ALL(:,:,1:end)); 
    
%% Some calculations of forcing in the different low, mid, high cloud combinations 
 switch load_type
        case 'cam6'
            
     otherwise
         
isave_plot_states_bar_chart=0;
isave_plot=0; %for UM_save_plot routine
savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
save_str = ['states_bar'];


    %Using PI clouds
    thresh_CF_states=0.01;
    thresh_CF_states=0.0001;    
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
    
    
 end 
    
    

    
%% 2D colour plots of mean contribution to forcing from different combinations of PI and PD
% cloud fraction using PI cloud fractions - start with no resctrictions to
% states? Or use combined states from PI and PD?
    clear gca
    isave_mean_2D_PDF=1;  %for paper draft
    
     switch load_type
        case 'cam6'
            
     otherwise
    
    
    
% Plot the forcing vs the cloud fraction in states 1 and 2 (1=clear -
% included to get zero CFs).

    %inds_PI_clear_low = unique( [inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    %inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    
    %I think that the way it is done below will allow non-clear/low points
    %to be included for one of the states
    %inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    
    %Whereas using intersect here should make it so that e.g. PI needs to
    %be in the clear state AND PD in the clear state (but using the 4
    %combinations :- clear + clear, clear + low, etc.)
    inds_CC = intersect( inds_PI{1}, inds_PD{1} );
    inds_CL = intersect( inds_PI{1}, inds_PD{2} );
    inds_LC = intersect( inds_PI{2}, inds_PD{1} );
    inds_LL = intersect( inds_PI{2}, inds_PD{2} );
    
    inds_PI_clear_low = unique( [inds_CC; inds_CL; inds_LC; inds_LL] );
    
    %run 2D histo template script :-    
    X_driver = cf_low(inds_PI_clear_low);
    Y_driver = cf_low_PD(inds_PI_clear_low);    
    Z_driver = forcing(inds_PI_clear_low);
    xlabelstr='Low cloud fraction Pre-Industrial';
    ylabelstr = 'Low cloud fraction Present Day';
    
%plot the mean forcing in each cloud fraction bin
    DRIVER_template_2D_mean_SW_forcing_vs_CF_UM  
    shading faceted; %this adds grid lines to the plot - but gives NaNs a
    %colour...
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]);
    label_text_pcolor; %Script to add text labels of the numbers for each block
    title('Mean forcing (W m^{-2})');
    save_str = [region_choice '_forcing_mean_2D_PDF_cf_PI_PD'];
    %savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir_date,[save_str] );  
        
%plot the overaell forcing (mean * frequency) in each cloud fraction bin    - i.e., these sum to the total forcing 
    %N.B. - the difference here is that the following scripts sets this
    %flag :- ioverall_contribution_mean_2D = 1;
    DRIVER_template_2D_overall_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis([-2 2]); 
    label_text_pcolor; %Script to add text labels of the numbers for each block
    %caxis([-2.5 2.5]);       
    title('Overall forcing (W m^{-2})');
    save_str = [region_choice '_forcing_overall_2D_PDF_cf_PI_PD']; 
    %savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir_date,[save_str] );  
    
    % As above, but in 'folded over' mode where add the net contrbutions
    % from opposite terms in the matrix
     
    %plot the overall forcing (mean * frequency) in each cloud fraction bin    - i.e., these sum to the total forcing 
    %N.B. - the difference here is that the following scripts sets this
    %flag :- ioverall_contribution_mean_2D = 1;
    ifold_over_ACSIS=1;
    %xlabelstr='Low cloud fraction Pre-Industrial';
    %ylabelstr = 'Low cloud fraction Present Day';
    DRIVER_template_2D_overall_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis([-2 2]); 
    label_text_pcolor; %Script to add text labels of the numbers for each block
    %caxis([-2.5 2.5]);       
    title('Net overall forcing (W m^{-2})');
    save_str = [region_choice '_NET_forcing_overall_2D_PDF_cf_PI_PD']; 
    %savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir_date,[save_str] );  
    
    
     end

%% 2D histograms of forcing vs cloud state    

 switch load_type
        case 'cam6'
            
     otherwise
    
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
    label_text_pcolor %add text labels
    save_str = [region_choice '_forcing_overall_2D_PDF_vs_cloud_states_PI_PD'];
    %savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir_date,[save_str] );  
    
  
    ifold_over_ACSIS=1;    
    DRIVER_template_2D_overall_SW_forcing_vs_cloudstate_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis([-2 2]);    
    %caxis([-2.5 2.5]);       
    title('Overall forcing (W m^{-2})');
    set(gca,'xlim',[0.5 8.5]);
    set(gca,'ylim',[0.5 8.5]);  
    label_text_pcolor %add text labels
    save_str = [region_choice '_NET_forcing_overall_2D_PDF_vs_cloud_states_PI_PD'];
    %savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir_date,[save_str] );  
    
    
 end
 
    
 
%% 2D histograms of cloud fraction vs Nd (Rosenfeld comparison)
clear gca
model_CTH_thresh = 3.2e3; 
%nmin=300; %min no. points for each bin
nmin=50; %min no. points for each bin

%CF type just for the y-axis
%CF_type = 'Modis COSP low CF';
CF_type = 'Model low CF';
%CF_type = 'Model tot CF';
%CF_type = 'Modis COSP liq CF';

%CF type just for the LWP binning
%CF_type_LWP = 'Modis COSP low CF';
CF_type_LWP = 'Model low CF';
%CF_type_LWP = 'Model tot CF';
%CF_type_LWP = 'Modis COSP liq CF';

restrict_height_method = 'max_cloud_height_in_cloud_LWC_IWC_PD';
%restrict_height_method = 'max_cloud_height_in_cloud_LWC_CFice_PD';
restrict_height_method = 'none';

bin_var = 'LWP';
bin_var = 'CTH';

bin_vals = 'Rosenfeld';
%bin_vals = 'Dan';
%bin_vals = 'Percentiles';


iocean_only=1;

iadd_conv_LWP=0; %Whether to include convective LWP
ilow_cloud_states_only=0; %whether to restrict to low-only model cloud scenes


var_UM = 'Land_mask'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
time_inds_landmask=1;
[Land_mask_PD,nT] = UM_get_time_data_mat_nc(dat_global,time_inds_landmask,load_type,gcm_Plat2D_UM); %1 is land, 0 is ocean
Land_mask_PD_ALL = repmat(Land_mask_PD,[1 1 length(time_inds)]); %replicate to the size of the 3D arrays

 switch load_type
        case 'cam6'
            
     otherwise
         
         
% included to get zero CFs).

    %inds_PI_clear_low = unique( [inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    %inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    
    %I think that the way it is done below will allow non-clear/low points
    %to be included for one of the states
    %inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    
    
    if ilow_cloud_states_only==1
        %Whereas using intersect here should make it so that e.g. PI needs to
        %be in the clear state AND PD in the clear state (but using the 4
        %combinations :- clear + clear, clear + low, etc.)
        inds_CC = intersect( inds_PI{1}, inds_PD{1} );
        inds_CL = intersect( inds_PI{1}, inds_PD{2} );
        inds_LC = intersect( inds_PI{2}, inds_PD{1} );
        inds_LL = intersect( inds_PI{2}, inds_PD{2} );
        
        inds_PI_clear_low = unique( [inds_CC; inds_CL; inds_LC; inds_LL] );
        
    else
        inds_PI_clear_low = [1:prod(size(LWP_PD_ALL))];
    end
    
    %run 2D histo template script :-       
        
    clear xvals_LWP
    clear yvals_LWP
    
    switch bin_var
        case 'LWP'            
            if iadd_conv_LWP==1
                LWP_Rosenfeld =  LWP_PD_ALL + Conv_LWP_PD_ALL + RWP_PD_ALL + Conv_RWP_PD_ALL;
            else
                LWP_Rosenfeld =  LWP_PD_ALL;
            end
            
            switch CF_type_LWP
                case 'Modis COSP low CF'
                    cf_Rosen_LWP = low_modisCF_PD_ALL;
                case 'Model low CF'
                    cf_Rosen_LWP = cf_low_PD;
                case 'Model tot CF'
                    cf_Rosen_LWP = cf_total_PD;
                case 'Modis COSP liq CF'
                    %cf_Rosen_LWP = liqCF_COSP_MODIS_PD_ALL; %MODIS liquid CF
                    cf_Rosen_LWP = liqCF_COSP_region;
            end
            
            LWP_in_cloud =  LWP_Rosenfeld ./ cf_Rosen_LWP;
            inan=find(cf_Rosen_LWP<0.02);
            LWP_in_cloud(inan)=NaN;
            
        case 'CTH'            
            LWP_in_cloud = max_cloud_height_in_cloud_LWC_IWC_PD/1e3;
            bin_vals = 'CTH';
            
    end
    
    switch CF_type
        case 'Modis COSP low CF'
            %cf_Rosen = low_modisCF_PD_ALL; 
            cf_Rosen = liqCF_COSP_region;
        case 'Model low CF'
            cf_Rosen = cf_low_PD; 
        case 'Model tot CF'            
            cf_Rosen = cf_total_PD;
        case 'Modis COSP liq CF'
            %cf_Rosen = liqCF_COSP_MODIS_PD_ALL; %MODIS liquid CF  
            cf_Rosen = liqCF_COSP_region;
    end
    
    switch restrict_height_method
        case 'max_cloud_height_in_cloud_LWC_IWC_PD'            
            inan = find(max_cloud_height_in_cloud_LWC_IWC_PD > model_CTH_thresh | isnan(max_cloud_height_in_cloud_LWC_IWC_PD) );
            LWP_in_cloud(inan)=NaN;
        case 'max_cloud_height_in_cloud_LWC_CFice_PD'
            inan=find(max_cloud_height_in_cloud_LWC_CFice_PD>model_CTH_thresh | isnan(max_cloud_height_in_cloud_LWC_CFice_PD) );
            LWP_in_cloud(inan)=NaN;
    end
    
    if iocean_only==1
        inan = find(Land_mask_PD_ALL==1);
        LWP_in_cloud(inan)=NaN;
    end
    
  
    
    
 
    switch bin_vals
        case 'Rosenfeld'
            LWP_bins=[0 75 120 170 240 300 400 600];
        case 'Dan'
            LWP_bins=[34 36];
            LWP_bins=[66 68];
            LWP_bins=[0 600];
        case 'Percentiles'
            lwp_max=6000; %Apply the same max LWP as in Rosenfeld
            nbins=7;
            dprc=100/nbins;
            prcs=[0 dprc:dprc:100];
            lwp_for_prc = LWP_in_cloud(inds_PI_clear_low);
            lwp_for_prc(lwp_for_prc>lwp_max)=NaN;
            LWP_bins = prctile(lwp_for_prc,prcs);
            
        case 'CTH'
            LWP_bins = [0:0.2:3.6];
            LWP_bins = [0:0.1:1.6];
            
    end    
        
    
    clear xvals_LWP yvals_LWP NY_vals_LWP
    for iLWP=1:length(LWP_bins)-1
        
        inan=find(LWP_in_cloud(inds_PI_clear_low)<=LWP_bins(iLWP) | LWP_in_cloud(inds_PI_clear_low)>LWP_bins(iLWP+1) | isnan(LWP_in_cloud(inds_PI_clear_low))==1 );
        
        X_driver = Nd_region(inds_PI_clear_low)/1e6;
        Y_driver = cf_Rosen(inds_PI_clear_low);        
        Y_driver(inan)=NaN;
        %Z_driver = forcing(inds_PI_clear_low);
        ylabelstr='Low cloud fraction';
        xlabelstr = 'Mean droplet Concentration (cm^{-3})';
        
        %plot the mean forcing in each cloud fraction bin
        DRIVER_template_2D_histo_CF_vs_Nd_UM
        shading faceted; %this adds grid lines to the plot - but gives NaNs a
        %colour...
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        %caxis([0 100]);
        %label_text_pcolor; %Script to add text labels of the numbers for each block
        title('CF vs CDNC histogram');
        
        close(gcf);
        
        xvals_LWP{iLWP} = mid_Xbins;
        yvals_LWP{iLWP} = Y_mean;
        NY_vals_LWP{iLWP} = NY_vals; %number of datapoints in each Nd bin
        
    end
    
    
%      cols='rbgkmrcymkgbrrrgykmc';
    cols_rgb{1}=[0 0 1];
    cols_rgb{2}=[1 0 1];
    cols_rgb{3}=[0 0.7 0.7];
    cols_rgb{4}=[1 1 0];
    cols_rgb{5}=[1 0.7 0];
    cols_rgb{6}=[1 0.3 0];
    cols_rgb{7}=[1 0 0];
    cols_rgb{8}=[0 0 1];
    cols_rgb{9}=[1 0 1];
    cols_rgb{10}=[0 0.7 0.7];
    cols_rgb{11}=[1 1 0];
    cols_rgb{12}=[1 0.7 0];
    cols_rgb{13}=[1 0.3 0];
    cols_rgb{14}=[1 0 0];  
    cols_rgb{15}=[0 0 1];
    cols_rgb{16}=[1 0 1];
    cols_rgb{17}=[0 0.7 0.7];
    cols_rgb{18}=[1 1 0];
    
   
    
    %patt={'-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--'};
    patt={'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
    
    figure
    set(gcf,'position',[20 210 600 384])    
     for iLWP=1:length(LWP_bins)-1
        xvals = xvals_LWP{iLWP};
        
        
        xvals(NY_vals_LWP{iLWP} < nmin) = NaN;
        
        %plot( xvals , yvals_LWP{iLWP} ,['x' cols(iLWP) patt{iLWP}],'linewidth',3);
        plot( xvals , yvals_LWP{iLWP} ,['o' patt{iLWP}],'linewidth',3,'color',cols_rgb{iLWP},'markerfacecolor',cols_rgb{iLWP});  
        hold on
        leg_str{iLWP}=['LWP=' num2str(LWP_bins(iLWP)) '-' num2str(LWP_bins(iLWP+1)) ];         
     end
     legend(leg_str,'location','northeastoutside');
     
     set(gca,'xscale','log');
     set(gca,'xlim',[10 300]);
     set(gca,'ylim',[0 1]);
     set(gca,'fontsize',16);
     xlabel(xlabelstr);
     ylabel(ylabelstr);
     titlenam = remove_character( ['UKCA (' um_case_PD '), ' CF_type ', ' CF_type_LWP, ', ' restrict_height_method ', iocean_only=' num2str(iocean_only)] , '_',' ');
     titwrapped = wrap_title_to_nlines(titlenam,50,3)
     title(titwrapped);
     grid on  
     
     save_str = [region_choice '_2D_PDF_cf_vs_Nd'];
     
     titlenam_driver = [um_case_PD '_Model_CFvsNd_' save_str];
     savename=[savedir_date titlenam_driver];
     clear opts
     %        opts.iplot_png=1;
     opts.iplot_eps=1;
     saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
    
    
    
    %savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    %UM_save_plot(gcf,isave_mean_2D_PDF,savedir_date,[save_str] );           
         
         
    
 end
 
 
 
%% Some map plots for the no-scavenging etc. runs with height restrction, etc.
iplot_noscav=0;
if iplot_noscav==1
Nd_PD_ALL2 = Nd_PD_ALL;
Nd_PI_ALL2 = Nd_PI_ALL;
Nd_PD_ALL2(inan)=NaN;
Nd_PI_ALL2(inan)=NaN;
Nd_PD_map2 = meanNoNan(Nd_PD_ALL2,3);
Nd_PI_map2 = meanNoNan(Nd_PI_ALL2,3);

dat_modis = (Nd_PD_map2 - Nd_PI_map2)/1e6; var_UM = 'Change in Nd (PD minus PI; cm^{-3})';
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-50 50]);


%COSP liquid CF
dat_PD_ALL2 = liqCF_COSP_MODIS_PD_ALL;
dat_PI_ALL2 = liqCF_COSP_MODIS_PI_ALL;
%dat_PD_ALL2(inan)=NaN;
%dat_PI_ALL2(inan)=NaN;
dat_PD_map2 = meanNoNan(dat_PD_ALL2,3);
dat_PI_map2 = meanNoNan(dat_PI_ALL2,3);

dat_modis = (dat_PD_map2 - dat_PI_map2); var_UM = 'Change in MODIS liquid CF (PD minus PI)';
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-0.1 0.1]);
    
end
 
%% 2D histogram for MODIS (Rosenfeld comparison)

    
    LWP_bins=[0 75 120 170 240 300 400 600];
    
    
    clear xvals_LWP
    clear yvals_LWP
    
    imodis_cf=0;

        LWP_in_cloud = modis_loaded_CF0.W_time3*1e3; %convert to g/m2
        %cf_modis = modis_loaded_CF0.Cloud_Fraction_Liquid.timeseries3;
        %cf_modis = modis_loaded_CF0.Cloud_Fraction_Combined.timeseries3; %Includes ice, so probably not useful.
        cf_modis = modis_loaded_CF0.Cloud_Fraction_Liquid.timeseries3 + modis_loaded_CF0.Cloud_Fraction_Undetermined.timeseries3; 
        inan=find(cf_modis<0.02);
        LWP_in_cloud(inan)=-1e9;
        
        [gcm_Plon2D_edges_AMSRE,gcm_Plat2D_edges_AMSRE] = get_edges_lat_lon(gcm_Plon2D_AMSRE,gcm_Plat2D_AMSRE);
        nT_modis = size(LWP_in_cloud,3);
        
[iregion_lin_modis,iregion_lin_edges_modis,cf_regional_modis]=get_lat_lon_irregular_with_time(nT_modis,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,gcm_Plat2D_edges_AMSRE,gcm_Plon2D_edges_AMSRE,cf_modis); 
[iregion_lin_modis,iregion_lin_edges_modis,Nd_regional_modis]=get_lat_lon_irregular_with_time(nT_modis,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,gcm_Plat2D_edges_AMSRE,gcm_Plon2D_edges_AMSRE,modis_loaded_CF0.Droplet_Number_Concentration_37.timeseries3);         
    
clear xvals_LWP yvals_LWP
    for iLWP=1:length(LWP_bins)-1
        
        inan=find(LWP_in_cloud<=LWP_bins(iLWP) | LWP_in_cloud>LWP_bins(iLWP+1) | isnan(LWP_in_cloud)==1);
        
        X_driver = Nd_regional_modis;
        Y_driver = cf_regional_modis;

        Y_driver(inan)=NaN;
        %Z_driver = forcing(inds_PI_clear_low);
        ylabelstr='Mean liquid low cloud fraction';
        xlabelstr = 'Droplet Concentration (cm^{-3})';
        
        %plot the mean forcing in each cloud fraction bin
        DRIVER_template_2D_histo_CF_vs_Nd_UM
        shading faceted; %this adds grid lines to the plot - but gives NaNs a
        %colour...
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        %caxis([0 100]);
        %label_text_pcolor; %Script to add text labels of the numbers for each block
        title('CF vs CDNC histogram');
        
        xvals_LWP{iLWP} = mid_Xbins;
        yvals_LWP{iLWP} = Y_mean;
        
    end
    
    figure
    cols='rbgkmyc';
    cols='rbgkmyc';
    cols_rgb{1}=[0 0 1];
    cols_rgb{2}=[1 0 1];
    cols_rgb{3}=[0 0.7 0.7];
    cols_rgb{4}=[1 1 0];
    cols_rgb{5}=[1 0.7 0];
    cols_rgb{6}=[1 0.3 0];
    cols_rgb{7}=[1 0 0];
    patt={'-','--','-','--','-','--','-','--'};
     for iLWP=1:length(LWP_bins)-1
        %plot( xvals_LWP{iLWP} , yvals_LWP{iLWP} ,['x' cols(iLWP) patt{iLWP}],'linewidth',3);
        plot( xvals_LWP{iLWP} , yvals_LWP{iLWP} ,['o' patt{iLWP}],'linewidth',3,'color',cols_rgb{iLWP},'markerfacecolor',cols_rgb{iLWP});        
        hold on
        leg_str{iLWP}=['LWP=' num2str(LWP_bins(iLWP)) '-' num2str(LWP_bins(iLWP+1)) ];         
     end
     legend(leg_str,'location','northeastoutside');
     
     set(gca,'xscale','log');
     set(gca,'xlim',[10 300]);
     set(gca,'ylim',[0 0.7]);
     set(gca,'fontsize',16);
     xlabel(xlabelstr);
     ylabel(ylabelstr);
     title('MODIS');
     grid on
     
    
    
    save_str = [region_choice '_2D_PDF_cf_vs_Nd'];
    %savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    
    titlenam_driver = ['MODIS_CFvsNd_' save_str];
     savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts); 
    


 
 
%% LWP comparison only in Sc gridboxes - NOT USING for the paper at the moment (used the ratio of LWP to LWP instead)
iplot_LWP_in_Sc=0;
if iplot_LWP_in_Sc==1
    
isave_plot_LWP_bias_Sc=1;

 switch load_type
        case 'cam6'
            
     otherwise

% Will first try just using the MODIS CTH<3.25km criteria for the AMSRE obs
% and the cloud states 1 and 2 (clear and low-only cloud) for the model
% Could also try restricting to CF>80%
% Prob can't expect model to replicate the real cloud fraction at every
% timestep, so try just doing separate time averages for model and obs
% Can only do in the daytime for MODIS if using the liquid cloud fraction
% Although the CTH should be usable for night too?
% Current CTH data perhaps only contains daytime data, though?
% Yes, because the mock L3 only contains SZA<65 deg data.
% So, should restrict the model to daytime too.

%AMSRE data

%amsre_LWP_time_period_mean = meanNoNan( meanNoNan(  amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,1),4) ,3);
%amsre_TLWP_time_period_mean = meanNoNan( meanNoNan(  amsre_tlwp(:,:,time_inds_AMSRE,1),4) ,3);     
   
 
    
%--- Run subplot script ------------------------------------------- 

    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only_DEFAULTS
            
    %Run for fraction Sc days
    ioverride_LWP_Sc_choices=1; 
    plot_var_Sc = 'Ndat';       
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
    
    %Run for LWP in Sc days
    ioverride_LWP_Sc_choices=1;       
    plot_var_Sc = 'LWP';
    UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
    
 end

    
 
end
%-----------------------------------------------------------------


% ----------------------------------------------------------------------------------
% ---- Partition the SW_down forcing (PI to PD) into changes due to CF, LWP and SW
% ----------------------------------------------------------------------------------



 
%% 2D histogram of LWP vs CF - MODIS - Probably should degrade to model resolution?
iplot_block=0;
if iplot_block==1
   
min_LWP=0.9e-3;
min_LWP=20e-3;

LWP_type = 'Grid-box mean';
%LWP_type = 'In-cloud';
iswap_CF_LWP=0; %Whether to plot CF vs LWP rather than LWP vs CF (swap X and Y).
  
amsre_lwp_type = 'TLWP';
%amsre_lwp_type = 'LWP';

modis_CTH_screening = '3.2km';
modis_CTH_screening = 'none';

ifilter_LWP_fraction=1; %Whether to filter out scenes with an LWP to LWP+RWP ratio lower than min_LWP_ratio = 0.99;
min_LWP_ratio = 0.99;
  

        %LWP data
        switch modis_CTH_screening
            case '3.2km';
                [CTH_amsre_CF0] = time_match_data(amsre_matlab_time(time_inds_AMSRE),date_num,modis_loaded_CF0.CTH.timeseries3);
            case 'none'
                [CTH_amsre_CF0] = time_match_data(amsre_matlab_time(time_inds_AMSRE),date_num,modis_loaded_CF0_2.CTH.timeseries3);
        end
        
        %Use the CF>0 data for screening for CTH>3.2km
        inan=find(isnan(CTH_amsre_CF0)==1);
        
        amsre_TLWP_CTH_CF0 = flipdim( amsre_tlwp(:,:,time_inds_AMSRE,1) , 1); %Flip the lat axis since AMSRE is oriented differently to MODIS      
        amsre_TLWP_CTH_CF0(inan)=NaN;
        
        amsre_LWP_CTH_CF0 = flipdim( amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,1) , 1); %Flip the lat axis since AMSRE is oriented differently to MODIS      
        amsre_LWP_CTH_CF0(inan)=NaN;
        
        
        lwp_ratio_amsre = amsre_LWP_CTH_CF0 ./ amsre_TLWP_CTH_CF0;
        lwp_ratio_amsre(amsre_TLWP_CTH_CF0<min_LWP)=1; %set ratio to one when have near clear-sky, so that it is included.
        inan2 = find(lwp_ratio_amsre < min_LWP_ratio);
        amsre_LWP_CTH_CF0(inan2)=NaN;
        amsre_TLWP_CTH_CF0(inan2)=NaN;
        
        % The following was a test to see why the TLWP values were lower
        % than the LWP ones - is because there are more NaNs for TLWP,
        % likely resulting from the function that converts the rain rate
        % etc to RWP
        %inan2=find( isnan(amsre_LWP_CTH_CF0)==1 | isnan(amsre_TLWP_CTH_CF0)==1 );
        %amsre_LWP_CTH_CF0(inan2)=NaN;        
        %amsre_TLWP_CTH_CF0(inan2)=NaN;
        
        switch modis_CTH_screening
            case '3.2km';
                %cf_modis = modis_loaded_CF0.Cloud_Fraction_Liquid.timeseries3;
                %cf_modis = modis_loaded_CF0.Cloud_Fraction_Combined.timeseries3; %Includes ice, so probably not useful.
                cf_modis = modis_loaded_CF0.Cloud_Fraction_Liquid.timeseries3;
            case 'none'
                cf_modis = modis_loaded_CF0_2.Cloud_Fraction_Liquid.timeseries3;
        end
        
        
        [cf_modis_amsre_times] = time_match_data(amsre_matlab_time(time_inds_AMSRE),date_num,cf_modis);

        
        switch LWP_type
            case 'Grid-box mean';
                %set a min LWP to allow a log LWP axis, but to keep the zero LWP bin on the
                %plot - but only if doing grid-box mean. For in-cloud we'll just
                %reject low CF points.
                
                amsre_TLWP_CTH_CF0(amsre_TLWP_CTH_CF0 < min_LWP) = min_LWP;
                amsre_LWP_CTH_CF0(amsre_LWP_CTH_CF0 < min_LWP) = min_LWP;
        end
        
        switch amsre_lwp_type
            case 'LWP'                
                LWP_in_cloud = amsre_LWP_CTH_CF0 ./ cf_modis_amsre_times;                                
                LWP_gridbox = amsre_LWP_CTH_CF0;
                
            case 'TLWP'
                LWP_in_cloud = amsre_TLWP_CTH_CF0 ./ cf_modis_amsre_times;                                
                LWP_gridbox = amsre_TLWP_CTH_CF0;               
                
        end
        
        inan=find(cf_modis_amsre_times<0.02);
        LWP_in_cloud(inan)=NaN;
        
        
        
        
        
        [gcm_Plon2D_edges_AMSRE,gcm_Plat2D_edges_AMSRE] = get_edges_lat_lon(gcm_Plon2D_AMSRE,gcm_Plat2D_AMSRE);
        nT_modis = size(cf_modis_amsre_times,3);

        %Restrict to region of interest
[iregion_lin_modis,iregion_lin_edges_modis,cf_regional_modis_amsre_times]=get_lat_lon_irregular_with_time(nT_modis,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,gcm_Plat2D_edges_AMSRE,gcm_Plon2D_edges_AMSRE,cf_modis_amsre_times); 
[iregion_lin_modis,iregion_lin_edges_modis,LWPic_regional_amsre]=get_lat_lon_irregular_with_time(nT_modis,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,gcm_Plat2D_edges_AMSRE,gcm_Plon2D_edges_AMSRE,LWP_in_cloud);         
[iregion_lin_modis,iregion_lin_edges_modis,LWP_regional_amsre]=get_lat_lon_irregular_with_time(nT_modis,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,gcm_Plat2D_edges_AMSRE,gcm_Plon2D_edges_AMSRE,LWP_gridbox);         




if iswap_CF_LWP==1
    iswap_xy=1; %Whether to plot CF vs LWP (swap X and Y).
end

        X_driver = cf_regional_modis_amsre_times;
        %Y_driver = LWPic_regional_amsre*1e3;
        switch LWP_type
            case 'Grid-box mean'
                Y_driver = LWP_regional_amsre*1e3;
                ylabelstr = 'Grid-box mean LWP from AMSR-E (g m^{-2})';
            case 'In-cloud'
                Y_driver = LWPic_regional_amsre*1e3;
                ylabelstr = 'In-cloud LWP from AMSR-E (g m^{-2})';
        end
        
        


        %Z_driver = forcing(inds_PI_clear_low);
        xlabelstr='Mean liquid low cloud fraction';       
        
        
        
        %plot the mean forcing in each cloud fraction bin
        DRIVER_template_2D_histo_LWP_vs_CF_UM
        %shading faceted; %this adds grid lines to the plot - but gives NaNs a
        %colour...
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        %caxis([0 100]);
        %label_text_pcolor; %Script to add text labels of the numbers for each block
        if iswap_CF_LWP==0
            set(gca,'yscale','log')
        end        
        title(['MODIS: LWP vs CF histogram for LAT=' num2str(LAT_val_DRIVER2(1)) ' to ' num2str(LAT_val_DRIVER2(2)) ' LON=' num2str(LON_val_DRIVER2(1)) ' to ' num2str(LON_val_DRIVER2(2))]);
                   
        
        LWP_vs_CF_Ymean_MODIS = Y_mean_accurate;
        LWP_vs_CF_Xmean_MODIS = X_mean_accurate;
        LWP_vs_CF_NX_MODIS = NX_vals;
        LWP_vs_CF_NY_MODIS = NY_vals;
        Xbins_LWP_vs_CF_mean_MODIS = mid_Xbins;
        Ybins_LWP_vs_CF_mean_MODIS = mid_Ybins;        
        
        
    



        

        
end
        
%% 2D histogram of LWP vs CF - model
iplot_block=0;
if iplot_block==1
    
clear gca
LWP_type = 'Grid-box mean';
%LWP_type = 'In-cloud';
iswap_CF_LWP=0; %Whether to plot CF vs LWP rather than LWP vs CF (swap X and Y).

min_LWP=0.9;
min_LWP=20;

min_LWP_ratio=0.99;

UM_cf_type = 'COSP-MODIS states COSP CF';
%UM_cf_type = 'COSP-MODIS states UM CF';
UM_cf_type = 'COSP-CALIPSO states COSP CF';
%UM_cf_type = 'COSP-CALIPSO states UM CF';
%UM_cf_type = 'UM states UM CF';
%UM_cf_type = 'UM states COSP CF';
%UM_cf_type = 'All states COSP CF';
%UM_cf_type = 'All states UM CF';
UM_cf_type = 'RWP ratio screening WITH convection, UM CF';

LWP_choice='Inc. conv LWP+RWP';
%LWP_choice='No conv LWP+RWP';

RWP_choice = 'Inc. RWP';
%RWP_choice = 'No RWP';

irestrict_to_day=1;

low_clear_CF_only_2d = NaN*ones(size(low_CF_PD_ALL));
mid_clear_CF_only_2d = NaN*ones(size(low_CF_PD_ALL));

switch UM_cf_type
    case 'COSP-MODIS states COSP CF'
        [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_modisCF_PD_ALL,mid_modisCF_PD_ALL,high_modisCF_PD_ALL,thresh_CF_states);
        inds_clear_low_ALL = unique( [inds_PD_ALL{1}; inds_PD_ALL{2};] );
        %inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_modisCF_PD_ALL(inds_clear_low_ALL); 
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_modisCF_PD_ALL(inds_clear_low_ALL); 
    case 'COSP-MODIS states UM CF'
        [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_modisCF_PD_ALL,mid_modisCF_PD_ALL,high_modisCF_PD_ALL,thresh_CF_states);
        inds_clear_low_ALL = unique( [inds_PD_ALL{1}; inds_PD_ALL{2};] );
        %inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_CF_PD_ALL(inds_clear_low_ALL); 
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_CF_PD_ALL(inds_clear_low_ALL);  
    case 'COSP-CALIPSO states COSP CF'
        [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_calipsoCF_PD_ALL,mid_calipsoCF_PD_ALL,high_calipsoCF_PD_ALL,thresh_CF_states);
        inds_clear_low_ALL = unique( [inds_PD_ALL{1}; inds_PD_ALL{2};] );
        %inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_calipsoCF_PD_ALL(inds_clear_low_ALL); 
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_calipsoCF_PD_ALL(inds_clear_low_ALL); 
    case 'COSP-CALIPSO states UM CF'
        [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_calipsoCF_PD_ALL,mid_calipsoCF_PD_ALL,high_calipsoCF_PD_ALL,thresh_CF_states);
        inds_clear_low_ALL = unique( [inds_PD_ALL{1}; inds_PD_ALL{2};] );
        %inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_CF_PD_ALL(inds_clear_low_ALL); 
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_CF_PD_ALL(inds_clear_low_ALL);
    case 'UM states UM CF'
        [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_CF_PD_ALL,mid_CF_PD_ALL,high_CF_PD_ALL,thresh_CF_states);
        inds_clear_low_ALL = unique( [inds_PD_ALL{1}; inds_PD_ALL{2};] );
        %inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_CF_PD_ALL(inds_clear_low_ALL); 
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_CF_PD_ALL(inds_clear_low_ALL); 
    case 'UM states COSP CF'
        [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_CF_PD_ALL,mid_CF_PD_ALL,high_CF_PD_ALL,thresh_CF_states);
        inds_clear_low_ALL = unique( [inds_PD_ALL{1}; inds_PD_ALL{2};] );
        %inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_calipsoCF_PD_ALL(inds_clear_low_ALL);  
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_calipsoCF_PD_ALL(inds_clear_low_ALL); 
    case 'All states UM CF'
        inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_CF_PD_ALL(inds_clear_low_ALL);
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_CF_PD_ALL(inds_clear_low_ALL); 
    case 'All states COSP CF'
        inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
        low_clear_CF_only_2d(inds_clear_low_ALL) = low_modisCF_PD_ALL(inds_clear_low_ALL);
        mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_modisCF_PD_ALL(inds_clear_low_ALL); 
    case 'RWP ratio screening WITH convection, UM CF'
         LWP_LS_CONV = LWP_PD_ALL + Conv_LWP_PD_ALL;
         TLWP_LS_CONV = LWP_LS_CONV + RWP_PD_ALL + Conv_RWP_PD_ALL;
         inan=find(TLWP_LS_CONV<0.1);
         TLWP_LS_CONV2 = TLWP_LS_CONV;
         TLWP_LS_CONV2(inan) = 1; %Include these points as they are (almost) clear-sky
         LWP_fraction = LWP_LS_CONV ./ TLWP_LS_CONV2;
         
         inds_clear_low_ALL = find(LWP_fraction>min_LWP_ratio); %i.e. mostly LWP
         low_clear_CF_only_2d(inds_clear_low_ALL) = low_CF_PD_ALL(inds_clear_low_ALL);
         mid_clear_CF_only_2d(inds_clear_low_ALL) = mid_CF_PD_ALL(inds_clear_low_ALL); 
         
         
    
end

if irestrict_to_day==1
    low_clear_CF_only_2d(inot_tol_LT)=NaN; %discard nighttime values and ones outside of local time tol set in script above
    mid_clear_CF_only_2d(inot_tol_LT)=NaN;
end


% Deal with LWP and TLWP
switch RWP_choice
    case 'Inc. RWP'
        TLWP_LS_ALL = LWP_PD_ALL + RWP_PD_ALL;
        TLWP_CONV_ALL = Conv_LWP_PD_ALL + Conv_RWP_PD_ALL;
    case 'No RWP'
        TLWP_LS_ALL = LWP_PD_ALL;
        TLWP_CONV_ALL = Conv_LWP_PD_ALL;        
end

switch LWP_choice
    case 'Inc. conv LWP+RWP'
        dat_lwp = TLWP_CONV_ALL + TLWP_LS_ALL;
    case 'No conv LWP+RWP'
        dat_lwp = TLWP_LS_ALL;
end

low_clear_LWP_only_2d = NaN*ones(size(dat_lwp)); %dat_lwp may include RWP and convective LWP
low_clear_LWP_only_2d(inds_clear_low_ALL) = dat_lwp(inds_clear_low_ALL);

if irestrict_to_day==1
    low_clear_LWP_only_2d(inot_tol_LT)=NaN; %discard nighttime values and ones outside of local time tol set in script above
end

switch LWP_type
    case 'Grid-box mean';
        %set a min LWP to allow a log LWP axis, but to keep the zero LWP bin on the
        %plot - but only if doing grid-box mean. For in-cloud we'll just
        %reject low CF points.        
        low_clear_LWP_only_2d(low_clear_LWP_only_2d < min_LWP) = min_LWP;
end
        
inan=find(low_clear_CF_only_2d<0.02);
%inan=find(low_clear_CF_only_2d<0.00002);
LWP_in_cloud =  low_clear_LWP_only_2d ./ low_clear_CF_only_2d;
LWP_in_cloud(inan)=NaN;

LWP_gridbox = low_clear_LWP_only_2d;
                
%Restrict to region of interest
[iregion_lin_modis,iregion_lin_edges_modis,cf_regional_UM]=get_lat_lon_irregular_with_time(nT_main,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,low_clear_CF_only_2d); 
[iregion_lin_modis,iregion_lin_edges_modis,cf_mid_regional_UM]=get_lat_lon_irregular_with_time(nT_main,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,mid_clear_CF_only_2d); 
[iregion_lin_modis,iregion_lin_edges_modis,LWPic_regional_UM]=get_lat_lon_irregular_with_time(nT_main,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,LWP_in_cloud);         
[iregion_lin_modis,iregion_lin_edges_modis,LWP_regional_UM]=get_lat_lon_irregular_with_time(nT_main,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,LWP_gridbox);         


if iswap_CF_LWP==1
    iswap_xy=1; %Whether to plot CF vs LWP (swap X and Y).
end

        X_driver = cf_regional_UM;
        %X_driver = cf_mid_regional_UM;
        
        %Y_driver = LWPic_regional_amsre*1e3;
        switch LWP_type
            case 'Grid-box mean'
                Y_driver = LWP_regional_UM;
                ylabelstr = 'Grid-box mean LWP from UM (g m^{-2})';
            case 'In-cloud'
                Y_driver = LWPic_regional_UM;
                ylabelstr = 'In-cloud LWP from UM (g m^{-2})';
        end
        
        


        %Z_driver = forcing(inds_PI_clear_low);
        xlabelstr='Mean liquid low cloud fraction';       
        
        
        
        %plot the mean forcing in each cloud fraction bin
        DRIVER_template_2D_histo_LWP_vs_CF_UM
        %shading faceted; %this adds grid lines to the plot - but gives NaNs a
        %colour...
        lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
        caxis([0 0.015]);
        %label_text_pcolor; %Script to add text labels of the numbers for each block
        if iswap_CF_LWP==0
            set(gca,'yscale','log')
        end
        title(['UM: LWP vs CF histogram for LAT=' num2str(LAT_val_DRIVER2(1)) ' to ' num2str(LAT_val_DRIVER2(2)) ' LON=' num2str(LON_val_DRIVER2(1)) ' to ' num2str(LON_val_DRIVER2(2))]);
        
        LWP_vs_CF_Ymean_UM = Y_mean_accurate;
        LWP_vs_CF_Xmean_UM = X_mean_accurate;
        LWP_vs_CF_NX_UM = NX_vals;
        LWP_vs_CF_NY_UM = NY_vals;
        Xbins_LWP_vs_CF_mean_UM = mid_Xbins;
        Ybins_LWP_vs_CF_mean_UM = mid_Ybins;

end
        
%% Mean model and MODIS LWP plot
iplot_block=0;
if iplot_block==1
    
figure
hold on
plot(Xbins_LWP_vs_CF_mean_UM,LWP_vs_CF_Ymean_MODIS,'ko-');
plot(Xbins_LWP_vs_CF_mean_UM,LWP_vs_CF_Ymean_UM,'-bs');
xlabel(xlabelstr);
ylabel(['Bin mean ' ylabelstr]);
clear legstr
legstr{1}='MODIS';        
legstr{2}='UM';
legend(legstr);
grid on

figure
hold on
biasUM = 100*(LWP_vs_CF_Ymean_UM./LWP_vs_CF_Ymean_MODIS - 1);
plot(Xbins_LWP_vs_CF_mean_UM,biasUM,'ko-');
xlabel(xlabelstr);
ylabel(['UM bias (%)']);
%clear legstr
%legstr{1}='MODIS';        
%legstr{2}='UM';
%legend(legstr);
grid on
end

%% 1D CF histogram comparison model vs MODIS
iplot_block=0;
if iplot_block==1
    
figure
hold on
plot(Xbins_LWP_vs_CF_mean_UM,LWP_vs_CF_NY_MODIS ./ sum(LWP_vs_CF_NY_MODIS),'ko-');
plot(Xbins_LWP_vs_CF_mean_UM,LWP_vs_CF_NY_UM/ sum(LWP_vs_CF_NY_UM),'-bs');
xlabel(xlabelstr);
ylabel('Normalized Frequency');
clear legstr
legstr{1}='MODIS';        
legstr{2}='UM';
legend(legstr);
grid on

end

%% 1D LWP histogram comparison model vs MODIS
iplot_block=0;
if iplot_block==1
    
figure
hold on
plot(Ybins_LWP_vs_CF_mean_UM,LWP_vs_CF_NX_MODIS ./ sum(LWP_vs_CF_NX_MODIS),'k-','linewidth',3);
plot(Ybins_LWP_vs_CF_mean_UM,LWP_vs_CF_NX_UM/ sum(LWP_vs_CF_NX_UM),'-b','linewidth',3);
xlabel(ylabelstr);
ylabel('Normalized Frequency');
clear legstr
legstr{1}='MODIS';        
legstr{2}='UM';
legend(legstr);
set(gca,'xscale','log');
set(gca,'yscale','log');
grid on
        

end

%% CF PDFs with different CF options (if saved)
iplot_block=0;
if iplot_block==1
    
if exist('NY_UM_ALL')
    figure
    plot(mid_Xbins,NY_UM,'bo-');
    hold on
    plot(mid_Xbins,NY_COSP,'rs-');
    plot(mid_Xbins,NY_UM_ALL,'ks--');
    plot(mid_Xbins,NY_UM_state_COSP_CF,'gp--');    
    
    
    clear leg_str
    leg_str{1} = 'COSP state, UM CF';
    leg_str{2} = 'COSP state, COSP CF';    
    leg_str{3} = 'UM state, UM CF';
    leg_str{4} = 'UM state, COSP CF';
    
    legend(leg_str);
        
end
end
%% CF PDFs with different CF options (if saved) - Using all grid-boxes and mean LWP
iplot_block=0;
if iplot_block==1
    
if exist('NY_UM_no_filter')
    figure
    plot(mid_Xbins,NY_UM_no_filter,'bo-');
    hold on
    plot(mid_Xbins,NY_COSP_no_filter,'rs-');
    %plot(mid_Xbins,NY_UM_ALL,'ks--');
    %plot(mid_Xbins,NY_UM_state_COSP_CF,'gp--');    
    grid on
    
    clear leg_str
    leg_str{1} = 'ALL gridboxes, UM CF';
    leg_str{2} = 'ALL gridboxes, COSP CF';        
    
    legend(leg_str);
        
end

end

%% Calcs of changes in cloud properties
if close_all_figs==1
    close all
    clear gca    
end



%for it_sw=1:1 %size(SW_down_PI_ALL,3)




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
i=find(isnan(N0)==1 & f0>=cf_min); W0(i)=NaN; %f0(i)=NaN;
i=find(isnan(W0)==1 & f0>=cf_min); N0(i)=NaN; %f0(i)=NaN;
i=find(isnan(N1)==1 & f1>=cf_min); W1(i)=NaN; %f1(i)=NaN;
i=find(isnan(W1)==1 & f1>=cf_min); N1(i)=NaN; %f1(i)=NaN;

%Special case where we have no CF in PI, but some cloud in PD - this would
%usually lead to the change in CF being ignored since the PD W and N would
%be NaN. So, set the PI W and N values in these cases to be equal to the
%PD values, so that this effect is incorporated into the CF change effect
i=find(f0<cf_min & f1>=cf_min); %If they are both <cf_min then they will both be treated as clear-sky in calc_SW_surf
W0_2 = W0; N0_2=N0;
W0_2(i) = W1(i); N0_2(i) = N1(i);

%Do a similar thing for when use PD as the baseline - set W and N to PI
%values if going from zero to some CF between PD and PI
i=find(f1<cf_min & f0>=cf_min);
W1_2 = W1; N1_2=N1;
W1_2(i) = W0(i); N1_2(i) = N0(i);


%% Calcs of changes in cloud properties and offline SW calcs

%Clear sky transmisison of the atmosphere (SW_down_surf / SW_down_TOA)
switch run_set
    case 'wind only nudging from level 18, reverted Hamish switches'
        transmission_atmos = 0.9; %Generally produces smaller magnitude (less negative) indirect forcing with smaller value
    otherwise
        transmission_atmos = 1; %model values suggest around 0.8
        transmission_atmos = 0.75; %model values suggest around 0.8       
end

%Can estimate this directly using the clean+clear radiation call
%(calculating in no aerosol situation as for the actual forcing)
SW_thresh=200;
transmission_atmos = SW_clean_clear_PD_ALL ./ SW_in_PD_ALL;
transmission_atmos(SW_in_PD_ALL<SW_thresh)=0.5;

SW_in_filter = SW_in_PD_ALL;
SW_surf_filter = SW_clean_clear_PD_ALL;
inan=find(SW_in_filter<SW_thresh);
SW_in_filter(inan)=NaN;
SW_surf_filter(inan)=NaN;

transmission_atmos2 =  SW_surf_filter ./ SW_in_filter;

transmission_atmos_mean = meanNoNan(SW_surf_filter,3) ./ meanNoNan(SW_in_filter,3);




dat_modis = meanNoNan(W1-W0,3); var_UM = 'Change in in-cloud LWP';
dat_modis = meanNoNan(100*(N1-N0)./N0,3); var_UM = '% Change in Nd';
it=1; dat_modis = meanNoNan(100*(N1(:,:,it)-N0(:,:,it))./N0(:,:,it),3); var_UM = '% Change in Nd';
it=1; dat_modis = meanNoNan(100*(W1(:,:,it)-W0(:,:,it))./W0(:,:,it),3); var_UM = '% Change in LWPic';
%it=1; dat_modis = meanNoNan(f1(:,:,it)-f0(:,:,it),3); var_UM = 'Change in CF';
%it=1; dat_modis = meanNoNan(SW_PD(:,:,it)-SW_cf2(:,:,it),3); var_UM = 'SW change due to CF';
%run plotting script
% figure
% UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
% lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([-20 10]);


SW_in = SW_in_PD_ALL(:,:,it_sw); %This is the TOA downward SW
    

%[Ac(i),tau(i),A_f(i),SW_f(i)] = calc_SW_surf(f1,W0,N0,SW_in,A_clear,transmission_atmos);
%[Ac(i),tau(i),A_W(i),SW_W(i)] = calc_SW_surf(f0,W1,N0,SW_in,A_clear,transmission_atmos);
%[Ac(i),tau(i),A_N(i),SW_N(i)] = calc_SW_surf(f0,W0,N1,SW_in,A_clear,transmission_atmos);

calc_set=2;
switch calc_set
    case 1
        [Ac,tau,T_f,SW_PI] = calc_SW_surf(f0,W0,N0,SW_in,transmission_atmos,cf_min);
        [Ac,tau,T_f,SW_PD] = calc_SW_surf(f1,W1,N1,SW_in,transmission_atmos,cf_min);
        
        [Ac,tau,T_f,SW_cf] = calc_SW_surf(f1,W0_2,N0_2,SW_in,transmission_atmos,cf_min);
        [Ac,tau,T_f,SW_lwp] = calc_SW_surf(f0,W1,N0,SW_in,transmission_atmos,cf_min);
        [Ac,tau,T_f,SW_Nd] = calc_SW_surf(f0,W0,N1,SW_in,transmission_atmos,cf_min);
        
        [Ac,tau,T_f,SW_cf2] = calc_SW_surf(f0,W1_2,N1_2,SW_in,transmission_atmos,cf_min);
        [Ac,tau,T_f,SW_lwp2] = calc_SW_surf(f1,W0,N1,SW_in,transmission_atmos,cf_min);
        [Ac,tau,T_f,SW_Nd2] = calc_SW_surf(f1,W1,N0,SW_in,transmission_atmos,cf_min);

case 2
    
    [Ac,tau,T_f,SW_PI] = calc_SW_surf(f0,W0_2,N0_2,SW_in,transmission_atmos,cf_min);
    [Ac,tau,T_f,SW_PD] = calc_SW_surf(f1,W1_2,N1_2,SW_in,transmission_atmos,cf_min);
    
    [Ac,tau,T_f,SW_cf] = calc_SW_surf(f1,W0_2,N0_2,SW_in,transmission_atmos,cf_min);
    [Ac,tau,T_f,SW_lwp] = calc_SW_surf(f0,W1_2,N0_2,SW_in,transmission_atmos,cf_min);
    [Ac,tau,T_f,SW_Nd] = calc_SW_surf(f0,W0_2,N1_2,SW_in,transmission_atmos,cf_min);
    
    [Ac,tau,T_f,SW_cf2] = calc_SW_surf(f0,W1_2,N1_2,SW_in,transmission_atmos,cf_min);
    [Ac,tau,T_f,SW_lwp2] = calc_SW_surf(f1,W0_2,N1_2,SW_in,transmission_atmos,cf_min);
    [Ac,tau,T_f,SW_Nd2] = calc_SW_surf(f1,W1_2,N0_2,SW_in,transmission_atmos,cf_min);

end


%% Get AMSRE in-cloud LWP for DAYTIME only - needed for the following calculations.
         ioverride_LWP_Sc_choices=1;
         plot_var_Sc = 'LWPic';
         lwp_day_night = 'day'; %Fairly weak dirunal cycle.
         iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
         irestrict_using_LWP_fraction=0; %whether to restrict to boxes with a high LWP to TLWP ratio
         %min_LWP_ratio=0.99;
         iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
         ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
         %mid or high cloudedit plo
         % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1
         cloud_scene_selection = 'none';
         irestrict_model_CF = 0;  %whether to restric the model CFs to > 80%
         imodis_cf = 0;
         modis_CTH_screening = 'none';
         LWP_sat = 'AMSR-E';
         fsize_latlon = 14;
         ibias_contour=1;
         bias_type = 'Percentage'; %overrides to absolute if plotting ratio of LWP to TLWP
         
         UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only
         
         lwpic_obs_mean_DAY = sat_data;         
         %also is thisw saved in the above script
         %lwpic_model_all_times = dat_modis_temp ./ low_calipsoCF_PD_ALL;
         % I.e. the in-cloud LWP for each time individually, rather than using
         % the time mean low CF to divide by
         lwpic_model_all_times_DAY_save = lwpic_model_all_times;
         
         
       
         
%% Do SW TOA calcs for the SW biases of model vs obs
isave_plot_SWTOA_calcs_draft=1;

diurnal_SW = 'average';
%diurnal_SW = 'day';

irestricted_LWP=0; 
%irestricted_LWP=1;

%First test that using the obs mean values (LWP_ic, Nd, CF) reproduces the observed SW TOA
% Also need the clear-sky albedo - can estimate from cloud free regions?
A_clear = 0.15; %Guess for now - value for CF quite sensitive to this (ranging up to 30%), so should prob
  %check what this is from the model
A_clear = 0.12; %value estimated from the 12 noon snapshot from VOCALS reginal run - N.B. - this uses the estimated transmission since
 %it is the surface albedo assuming no atmosphere above
%A_clear = 0.0;  
 %A_clear = 0.09;
 
 fac_tr = 1/1.1;
 fac_tr = 1;
 
 fac_SWin=1.7;
 fac_SWin=1;
 
 fac_Nd = 1;
 %fac_Nd = 2;%1.3;
 
 fac_cf = 1;
 %fac_cf = 1.2;
 
 fac_SW_calc=1;
 %fac_SW_calc=1.4;
 
 switch diurnal_SW
     case 'average'
         %re-grid all to the model resolution for convienience                  
         
         % Low CALIPSO CF
         dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO
         % Total CALIPSO CF
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO
                  
         %lwpic_obs_mean = 1e3 * griddata(modis_loaded.gcm_Plat2D_AMSRE,modis_loaded.gcm_Plon2D_AMSRE,amsre_TLWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
         % This calculated in the in-cloud LWP section
         
         if irestricted_LWP==0
             lwpic_obs_mean2 = lwpic_obs_mean;
         else
             lwpic_obs_mean2 = lwpic_obs_mean_restrict;
         end
         

         
     case 'day'
         % Low CALIPSO CF
         dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO 
         % Total CALIPSO CF
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO
         
            %dat_modis used for the CF average below
         lwpic_obs_mean2 = lwpic_obs_mean_DAY;
         
         if irestricted_LWP==0
             lwpic_obs_mean2 = lwpic_obs_DAY_mean;
         else
             lwpic_obs_mean2 = lwpic_obs_mean_DAY_restrict; %doesn't exist at present
         end
         
 end
 
 cf_obs_mean = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
 
 %modis_loaded is the CF.80 data
 %Nd_obs_mean = fac_Nd .* griddata(modis_loaded.gcm_Plat2D_AMSRE,modis_loaded.gcm_Plon2D_AMSRE,Nd37_1deg_time_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);  
 %Data saved from UM_ACSIS_global_SUBPLOT_commands_LWP_Sc_only.m 
 % N.B - has the trend correction applied :-
 Nd_obs_mean = fac_Nd .* Nd_sat_data_saved;
 
 tr_atmos_mean = meanNoNan(transmission_atmos(:,:,it_sw),3)* fac_tr;
 tr_atmos_mean = transmission_atmos_mean;
 %tr_atmos_mean = 0.8;
 SW_in_CTH_mean = fac_SWin * meanNoNan(SW_in .* transmission_atmos(:,:,it_sw),3); %N.B. - this should be the CLOUD TOP SW down for the TOA calcluation (not TOA).
 %SW_in_CTH_mean = fac_SWin * meanNoNan(SW_in ,3);
 
 %Using CERES SW in TOA
 SW_in_ceres_regrid = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,ceres_SW_TOA_in_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM); 
 SW_in_CTH_mean = fac_SWin .* SW_in_ceres_regrid .* tr_atmos_mean;
 
 % So, multiply by the transmission
 % Also, should get the CERES values for this....
[Ac_calc_obs,tau_calc_obs,A_calc_obs,SWTOA_calc_obs] = calc_SW(cf_obs_mean,lwpic_obs_mean2,Nd_obs_mean,SW_in_CTH_mean,A_clear,tr_atmos_mean,cf_min); 
SWTOA_calc_obs  = SWTOA_calc_obs * fac_SW_calc;
SWTOA_calc_obs(isnan(lwpic_obs_mean2)==1) = NaN; %mask out the land regions where there are no AMSRE retrievals

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA calculated from obs mean values using ' diurnal_SW ' CF';]
figure
dat_modis = SWTOA_calc_obs;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end

% irestrict_domain_DRIVER=0;
% var_UM = 'CF obs mean values';
% figure
% dat_modis = cf_obs_mean;
% UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
% lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));



SW_TOA_ceres_mean = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,ceres_SW_TOA_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM); 
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA percentage bias of calc';
figure
dat_modis = 100* (SWTOA_calc_obs ./ SW_TOA_ceres_mean - 1);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-50 50]);

um_data = SWTOA_calc_obs;
sat_data = SW_TOA_ceres_mean;
var_Latex='SWceres';
DRIVER_calc_biases_for_regions




%%
ceres_SW_TOA_in_time_period_mean = meanNoNan(SW_TOA_in_ceres(time_inds_CERES,:,:),1);
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA in CERES from obs mean values';
figure
dat_modis = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,ceres_SW_TOA_in_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM); 
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([200 450]);




irestrict_domain_DRIVER=0;
var_UM = 'LWP obs mean values';
figure
dat_modis = lwpic_obs_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));

irestrict_domain_DRIVER=0;
var_UM = 'Nd obs mean values';
figure
dat_modis = Nd_obs_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));

irestrict_domain_DRIVER=1;
var_UM = 'Transmission model mean values';
figure
dat_modis = tr_atmos_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0.7 0.85]);

irestrict_domain_DRIVER=1;
var_UM = 'SW in cloud top model mean values';
figure
dat_modis = SW_in_CTH_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));

irestrict_domain_DRIVER=1;
var_UM = 'SW in model mean values';
figure
dat_modis = meanNoNan(SW_in,3);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([200 450]);
  

%% Trying using CERES monthly CF and optical depth

tr_ceres_mean = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,transmission_atmos_mean,gcm_Plat2D_CERES,gcm_Plon2D_CERES); 
tr_ceres_rep = repmat(tr_ceres_mean,[1 1 size(cf_ceres(time_inds_CERES,:,:),1)]);
tr_ceres_rep = permute(tr_ceres_rep,[3 1 2]);
%tr_ceres_rep = 0.75;

f_tau = 1;

[Ac_calc_ceres,tau_calc_ceres,A_calc_obs,SWTOA_calc_ceres] = calc_SW(cf_ceres(time_inds_CERES,:,:)/100,[],[],SW_TOA_in_ceres(time_inds_CERES,:,:),A_clear,tr_ceres_rep,cf_min,1,f_tau.*tau_ceres(time_inds_CERES,:,:)); 

SWTOA_calc_ceres_mean = meanNoNan(SWTOA_calc_ceres,1);
SWTOA_calc_ceres_mean2 = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,SWTOA_calc_ceres_mean,gcm_Plat2D_UM,gcm_Plon2D_UM); 

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA ceres monthly calc';
figure
dat_modis = SWTOA_calc_ceres_mean2;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);


ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA percentage bias of calc using CERES monthly means';
figure
dat_modis = 100* ( SWTOA_calc_ceres_mean2 ./ SW_TOA_ceres_mean - 1);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);

it_plot=1;
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA ceres monthly calc SINGLE TIME';
figure
dat_modis = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,squeeze(SWTOA_calc_ceres(it_plot,:,:)),gcm_Plat2D_UM,gcm_Plon2D_UM); 
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

it_plot=1;
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'CERES CF, single time';
figure
dat_modis = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,squeeze(cf_ceres(time_inds_CERES(it_plot),:,:))/100,gcm_Plat2D_UM,gcm_Plon2D_UM); 
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 1]);

it_plot=1;
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'CERES tau, single time';
figure
dat_modis = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,squeeze(tau_ceres(time_inds_CERES(it_plot),:,:)),gcm_Plat2D_UM,gcm_Plon2D_UM); 
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([0 200]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'CERES tau, all times';
figure
%tau = meanNoNan(100*tau_ceres(time_inds_CERES,:,:) ./ cf_ceres(time_inds_CERES,:,:),1);
tau = meanNoNan(tau_ceres(time_inds_CERES,:,:),1);
dat_modis = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,tau,gcm_Plat2D_UM,gcm_Plon2D_UM); 
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([0 200]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'Calculated tau from CALIPSO, MODIS and AMSR';
figure
dat_modis = tau_calc_obs;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([0 200]);



it_plot=1;
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'CALIPSO low cloud';
figure
dat_modis = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,squeeze(cllcalipso_monthly_AVERAGE(time_inds_CAL(it_plot),:,:))/100,gcm_Plat2D_UM,gcm_Plon2D_UM); 
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 1]);

it_plot=1;
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'CALIPSO total cloud';
figure
dat_modis = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,squeeze(cltcalipso_monthly_AVERAGE(time_inds_CAL(it_plot),:,:))/100,gcm_Plat2D_UM,gcm_Plon2D_UM); 
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 1]);

%% Using time averages of the CERES data for a test

cf_ceres_mean = meanNoNan(cf_ceres(time_inds_CERES,:,:),1);
sw_in_ceres_mean = meanNoNan(SW_TOA_in_ceres(time_inds_CERES,:,:),1);
tau_ceres_mean = meanNoNan(tau_ceres(time_inds_CERES,:,:),1);


[Ac_calc_obs,tau_calc_obs,A_calc_obs,SWTOA_calc_ceres_using_means] = calc_SW(cf_ceres_mean/100,[],[],sw_in_ceres_mean,A_clear,tr_ceres_mean,cf_min,1,tau_ceres_mean); 



SWTOA_calc_ceres_using_means_regrid = griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,SWTOA_calc_ceres_using_means,gcm_Plat2D_UM,gcm_Plon2D_UM); 

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA ceres calc using means';
figure
dat_modis = SWTOA_calc_ceres_using_means_regrid;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);


ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA percentage bias of calc using CERES monthly means';
figure
dat_modis = 100* ( SWTOA_calc_ceres_using_means_regrid ./ SW_TOA_ceres_mean - 1);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);


%% Using model CF, LWP and Nd to calc TOA - options for using all data, or the mean
%% CF change calculation done here too.
time_mean=0;
time_mean=1;

nd_type='all times';
nd_type='CTH restriction';

%lwpic_model_all_times; %LWPic where used the instantaneous CF to divide by
%lwpic_model - used the time mean CF to divide by


switch time_mean
    case 0        
        lwp =  lwpic_model_all_times_restrict_save;
        %lwp = lwpic_model_all_times;        
        cf = low_calipsoCF_PD_ALL;        
        nd = Nd_PD_ALL;        
        %sw = SW_in;
        
        str_type_calc = 'instantaneous';
        trans = transmission_atmos2;
        trans = 0.75;
        
        sw = fac_SWin * SW_in .* trans;

        
    case 1
        if irestricted_LWP==0
            lwp = 1.0* lwpic_model;
        else
            lwp = 1.0* lwpic_model_restrict;
        end
                      
        cf = 1.0* meanNoNan(low_calipsoCF_PD_ALL,3);  
        switch nd_type
            case 'all times'
                nd = 1.0* meanNoNan(Nd_PD_ALL,3);
            case 'CTH restriction'
                nd = 1.0* meanNoNan(Nd_CTH_3pt2km_CF80,3)*1e6;
        end
             
        trans = transmission_atmos_mean;           
        str_type_calc = 'time-mean';                        
        
        %sw = meanNoNan(SW_in,3); 
        sw = fac_SWin * meanNoNan(SW_in ,3) .* trans;
        
end



[Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model] = calc_SW(cf,lwp,nd/1e6,sw,A_clear,trans,cf_min); 

SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,3);



ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA model ' str_type_calc ' calc'];
figure
dat_modis = SWTOA_calc_model_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA percentage bias of calc using ' str_type_calc ' model fields vs actual model SW TOA'];
figure
dat_modis = 100* ( SWTOA_calc_model_mean ./ SW_TOA_model_save - 1);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);

um_data = SWTOA_calc_model_mean;
sat_data = SW_TOA_model_save;
var_Latex='SWmodel';
DRIVER_calc_biases_for_regions

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['Transmissivity'];
figure
dat_modis = trans;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0.6 1]);


% Now calculate the SW effect of switching in the observed CF, LWPic and Nd in one-at-a-time tests
%Just doing using the time-mean values
diurnal_SW='day';
diurnal_SW='average';

nd_type = 'all times';
nd_type = 'CTH restriction';

%Stuff that won't change
trans = transmission_atmos_mean;           
sw = fac_SWin * meanNoNan(SW_in ,3) .* trans;
%sw = meanNoNan(SW_in,3);

str_type_calc = 'time-mean';

% ------ Change cf ------
% Set lwp and nd as model values
%cf = 1.0* meanNoNan(low_calipsoCF_PD_ALL,3);
if irestricted_LWP==0
    lwp = 1.0* lwpic_model;
    restrict_str = '';
else
    lwp = 1.0* lwpic_model_restrict;
    restrict_str = 'restricted LWPic';
end
switch nd_type
    case 'all times'
        nd = 1.0* meanNoNan(Nd_PD_ALL,3);
    case 'CTH restriction'
        nd = 1.0* meanNoNan(Nd_CTH_3pt2km_CF80,3)*1e6;
end


%Set to observed cf
switch diurnal_SW
     case 'average'
         %re-grid all to the model resolution for convienience
         dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO

         %lwpic_obs_mean2 = lwpic_obs_mean;
         
     case 'day'
         dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO 
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPS
            %dat_modis used for the CF average below
         %lwpic_obs_mean2 = lwpic_obs_mean_DAY;
         
 end
 
 cf = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
 %Nd_obs_mean = fac_Nd .* griddata(modis_loaded.gcm_Plat2D_AMSRE,modis_loaded.gcm_Plon2D_AMSRE,Nd37_1deg_time_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);  
 
%do the calc
[Ac_calc_model_cf,tau_calc_model_cf,A_calc_model_cf,SWTOA_calc_model_cf] = calc_SW(cf,lwp,nd/1e6,sw,A_clear,trans,cf_min); 
SWTOA_calc_model_mean_cf = meanNoNan(SWTOA_calc_model_cf,3);

%plot absolute
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to CF bias (W m^{-2})'];
figure
dSWTOA_calc_model_mean_cf = SWTOA_calc_model_mean - SWTOA_calc_model_mean_cf;
dat_modis = dSWTOA_calc_model_mean_cf;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);

%plot percentage
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to CF bias (%) ' restrict_str];
figure
%dat_modis = 100*(SWTOA_calc_model_mean./SWTOA_calc_model_mean_cf - 1);
dat_modis = 100*(dSWTOA_calc_model_mean_cf./ SW_TOA_ceres_mean); %expressing as a percentage of the CERES (true) value

UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end
    

%% Extra plots for CF TOA

%plot SW fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA calc with cf change (W m^{-2})'];
figure
dat_modis = SWTOA_calc_model_mean_cf;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

%plot SW fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA calculated from model fields (W m^{-2})'];
figure
dat_modis = SWTOA_calc_model_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

%% ------ Change LWPic ------

% Set lwp and nd as model values
cf = 1.0* meanNoNan(low_calipsoCF_PD_ALL,3);
%lwp = 1.0* lwpic_model;
switch nd_type
    case 'all times'
        nd = 1.0* meanNoNan(Nd_PD_ALL,3);
    case 'CTH restriction'
        nd = 1.0* meanNoNan(Nd_CTH_3pt2km_CF80,3)*1e6;
end

nd_model = nd;

%Set to observed cf
switch diurnal_SW
     case 'average'
         %re-grid all to the model resolution for convienience
         %dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO

         if irestricted_LWP==0
             lwp = lwpic_obs_mean;
             restrict_str = '';
         else
             lwp = lwpic_obs_mean_restrict;
             restrict_str = 'restricted LWPic';
         end
         
     case 'day'
         %dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO 
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPS
            %dat_modis used for the CF average below
         lwp = lwpic_obs_mean_DAY;
         
 end
 
 %cf = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
 %nd = fac_Nd .* griddata(modis_loaded.gcm_Plat2D_AMSRE,modis_loaded.gcm_Plon2D_AMSRE,Nd37_1deg_time_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);  
 
%do the calc
[Ac_calc_model_lwp,tau_calc_model_lwp,A_calc_model_lwp,SWTOA_calc_model_lwp] = calc_SW(cf,lwp,nd/1e6,sw,A_clear,trans,cf_min); 
SWTOA_calc_model_mean_lwp = meanNoNan(SWTOA_calc_model_lwp,3);
SWTOA_calc_model_mean_lwp(isnan(lwp)==1)=NaN;


%plot absolute
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to LWPic bias (W m^{-2})'];
figure
dSWTOA_calc_model_mean_lwp = SWTOA_calc_model_mean - SWTOA_calc_model_mean_lwp;
dat_modis = dSWTOA_calc_model_mean_lwp;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);

%plot percentage
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to LWPic bias (%) ' restrict_str];
figure
%dat_modis = 100*(SWTOA_calc_model_mean./SWTOA_calc_model_mean_lwp - 1);
dat_modis = 100*(dSWTOA_calc_model_mean_lwp./ SW_TOA_ceres_mean); %expressing as a percentage of the CERES (true) value
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end


%% Extra LWP plots

%plot SW fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA calc with LWPic change (W m^{-2})'];
figure
dat_modis = SWTOA_calc_model_mean_lwp;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

%plot LWP fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['LWP obs (g m^{-2})'];
figure
dat_modis = lwpic_obs_mean;
%dat_modis = lwp;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 1000]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['LWP model (g m^{-2})'];
figure
dat_modis = lwpic_model;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 1000]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['LWP model minus obs (g m^{-2})'];
figure
dat_modis = lwpic_model - lwpic_obs_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-200 200]);


%% ------ Change Nd ------

% Set lwp and cf as model values
cf = 1.0* meanNoNan(low_calipsoCF_PD_ALL,3); %MODEL COSP CALIPSO CF

%Model LWP:-
if irestricted_LWP==0
    lwp = 1.0* lwpic_model;
    restrict_str = '';
else
    lwp = 1.0* lwpic_model_restrict;
    restrict_str = 'restricted LWPic';
end

% switch nd_type
%     case 'all times'
%         nd = 1.0* meanNoNan(Nd_PD_ALL,3);
%     case 'CTH restriction'
%         nd = 1.0* meanNoNan(Nd_CTH_3pt2km_CF80,3)*1e6;
% end

% OBSERVED Nd :-

%nd = fac_Nd .* 1e6*griddata(modis_loaded.gcm_Plat2D_AMSRE,modis_loaded.gcm_Plon2D_AMSRE,Nd37_1deg_time_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);  
nd = 1e6*Nd_obs_mean; %Use the same obs Nd as for the overall SW calc (includes the trend correction)
%do the calc
[Ac_calc_model_Nd,tau_calc_model_Nd,A_calc_model_Nd,SWTOA_calc_model_Nd] = calc_SW(cf,lwp,nd/1e6,sw,A_clear,trans,cf_min); 
SWTOA_calc_model_mean_Nd = meanNoNan(SWTOA_calc_model_Nd,3);

%plot absolute
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to Nd bias (W m^{-2})'];
figure
dSWTOA_calc_model_mean_Nd = SWTOA_calc_model_mean - SWTOA_calc_model_mean_Nd;
dat_modis = dSWTOA_calc_model_mean_Nd;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-15 15]);

%plot percentage
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to Nd bias (%) ' restrict_str];
figure
%dat_modis = 100*(SWTOA_calc_model_mean./SWTOA_calc_model_mean_Nd - 1);
dat_modis = 100*(dSWTOA_calc_model_mean_Nd./ SW_TOA_ceres_mean); %expressing as a percentage of the CERES (true) value
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-5 5]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end


%% Extra plots for checking the Nd calcs

%plot SW fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA calc with Nd change (W m^{-2})'];
figure
dat_modis = SWTOA_calc_model_mean_Nd;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);


%plot Nd fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['Nd obs (cm^{-3})'];
figure
dat_modis = nd/1e6;
%dat_modis = lwp;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 300]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['Nd model (cm^{-3})'];
figure
dat_modis = nd_model/1e6;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 300]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['Nd model minus obs (cm^{-3})'];
figure
dat_modis = nd_model/1e6 - nd/1e6;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-200 200]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['% Nd vs obs (cm^{-3})'];
figure
dat_modis = 100*(nd_model./nd - 1);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

%% ------ Combine the contributions ------
dSWTOA_tot_cont = dSWTOA_calc_model_mean_cf + dSWTOA_calc_model_mean_lwp + dSWTOA_calc_model_mean_Nd;

%plot absolute
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to linear combination (W m^{-2})'];
figure
dat_modis = dSWTOA_tot_cont;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);



%plot percentage
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to linear combination (%)'];
figure
%dat_modis = 100*( dSWTOA_tot_cont ./ SWTOA_calc_model_mean );
dat_modis = 100*(dSWTOA_tot_cont ./ SW_TOA_ceres_mean); %expressing as a percentage of the CERES (true) value
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-50 50]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end



%% Change all 3 to obs for testing
% Set lwp and nd as model values
%cf = 1.0* meanNoNan(low_calipsoCF_PD_ALL,3);
%lwp = 1.0* lwpic_model;
%nd = 1.0* meanNoNan(Nd_PD_ALL,3); %per m3

%Set to observed cf
switch diurnal_SW
     case 'average'
         %re-grid all to the model resolution for convienience
         dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO       
         
         if irestricted_LWP==0
             lwp = lwpic_obs_mean;
         else
             lwp = lwpic_obs_mean_restrict;
         end
         
     case 'day'
         dat_modis = fac_cf.*meanNoNan(cllcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPSO 
         %dat_modis = fac_cf.*meanNoNan(cltcalipso_monthly_DAYTIME2(time_inds_CAL,:,:),1)/100; %time mean low CF from CALIPS
            %dat_modis used for the CF average below
         lwp = lwpic_obs_mean_DAY;
         
 end
 
 cf = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
 nd = fac_Nd .* 1e6*griddata(modis_loaded.gcm_Plat2D_AMSRE,modis_loaded.gcm_Plon2D_AMSRE,Nd37_1deg_time_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);  
 
%do the calc
[Ac_calc_model_lwp,tau_calc_model_lwp,A_calc_model_lwp,SWTOA_calc_model_all3] = calc_SW(cf,lwp,nd/1e6,sw,A_clear,trans,cf_min); 
SWTOA_calc_model_mean_all3 = meanNoNan(SWTOA_calc_model_all3,3);

%plot absolute
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to calc with all 3 biases (W m^{-2})'];
figure
dSWTOA_calc_model_mean_all3 = SWTOA_calc_model_mean - SWTOA_calc_model_mean_all3;
dat_modis = dSWTOA_calc_model_mean_all3;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);

%plot percentage
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA perturbation due to calc with all 3 biases (%)'];
figure
%dat_modis = 100*(SWTOA_calc_model_mean./SWTOA_calc_model_mean_lwp - 1);
dat_modis = 100*(dSWTOA_calc_model_mean_all3./ SW_TOA_ceres_mean); %expressing as a percentage of the CERES (true) value
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-30 30]);


%plot SW fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA CERES observed (W m^{-2})'];
figure
dat_modis = SW_TOA_ceres_mean;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

%plot SW fields for testing
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['SW TOA calc using all 3 observed fields (W m^{-2})'];
figure
dat_modis = SWTOA_calc_model_mean_all3;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);


%% abs contributions
%plot percentage conribution to sum of abs changes from each term
dSWTOA_tot_abs_cont = abs(dSWTOA_calc_model_mean_cf) + abs(dSWTOA_calc_model_mean_lwp) + abs(dSWTOA_calc_model_mean_Nd);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['Contribution to abs SW TOA bias from CF bias (%)'];
figure
%dat_modis = 100*( dSWTOA_tot_cont ./ SWTOA_calc_model_mean );
dat_modis = 100*(dSWTOA_calc_model_mean_cf ./ dSWTOA_tot_abs_cont); %expressing as a percentage of the CERES (true) value
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end

%savename = [savedir_date 'abs_SW_cont_bias_calc_CF'];
um_data = dSWTOA_calc_model_mean_cf;
sat_data = dSWTOA_tot_abs_cont;
var_Latex='SWabsCF';
DRIVER_calc_biases_for_regions


ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['Contribution to abs SW TOA bias from LWPic bias (%)'];
figure
%dat_modis = 100*( dSWTOA_tot_cont ./ SWTOA_calc_model_mean );
dat_modis = 100*(dSWTOA_calc_model_mean_lwp ./ dSWTOA_tot_abs_cont); %expressing as a percentage of the CERES (true) value
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end

%savename = [savedir_date 'abs_SW_cont_bias_calc_lwp'];
um_data = dSWTOA_calc_model_mean_lwp;
sat_data = dSWTOA_tot_abs_cont;
var_Latex='SWabsLWP';
DRIVER_calc_biases_for_regions

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = ['Contribution to abs SW TOA bias from Nd bias (%)'];
figure
%dat_modis = 100*( dSWTOA_tot_cont ./ SWTOA_calc_model_mean );
dat_modis = 100*(dSWTOA_calc_model_mean_Nd ./ dSWTOA_tot_abs_cont); %expressing as a percentage of the CERES (true) value
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);
if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end

%savename = [savedir_date 'abs_SW_cont_bias_calc_Nd'];
um_data = dSWTOA_calc_model_mean_Nd;
sat_data = dSWTOA_tot_abs_cont;
var_Latex='SWabsNd';
DRIVER_calc_biases_for_regions

%%
ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'SW TOA incoming from model';
figure
dat_modis = meanNoNan(SW_in,3);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 800]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'LWPic model from CF at each time step';
figure
dat_modis = meanNoNan(lwpic_model_all_times_save,3);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'LWPic model using time-av CF';
figure
dat_modis = lwpic_model;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);


ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'LWPic restricted model from CF at each time step';
figure
dat_modis = meanNoNan(lwpic_model_all_times_restrict_save,3);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 350]);

ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'LWPic restricted model using time-av CF';
figure
dat_modis = lwpic_model_restrict;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 350]);


ioverride_LAT_plots=0;
irestrict_domain_DRIVER=1;
var_UM = 'Tau calc model';
figure
dat_modis = meanNoNan(tau_calc_model,3);
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);



%% PI to PD difference plots for cf, LWP and Nd
isave_plot_global_LWP=1; %for draft paper

if close_all_figs==1
    close all
    clear gca    
end


    
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
%         saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
% %        close(gcf);
%     end 
    
    %dat_modis = 100 * (meanNoNan(f1,3)./meanNoNan(f0,3) - 1); var_UM = 'Percentage change in cloud fraction (PD minus PI)';    
    dlowcf_ALL = f1-f0;
    %dat_modis = meanNoNan(f1,3) - meanNoNan(f0,3); var_UM = 'Absolute change in cloud fraction (PD minus PI)';     
    dat_modis = 100*(meanNoNan(f1,3) - meanNoNan(f0,3)) ./ meanNoNan(f0,3); var_UM = 'Percentage change in cloud fraction (PD minus PI)';         
    dlowcf_timemean=dat_modis;
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-0.1 0.1]);   
    caxis([-30 30]);      
    
      % Save
    if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date titlenam_driver];        
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
    
    if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date titlenam_driver];        
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
    caxis([-150 150]);   
    
      % Save      
   if isave_plot_global_LWP==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
   end 




%% Do plots of the above paritioning results - not needed for 2d histos
isave_plots_partioning=1;  %used in draft paper
isave_plots_partioning2=0; %not used in draft at the moment
clear gca

isave_plot=0;

var_UM = 'Total estimated change in surface SW';
SW_estimated_forcing_total = meanNoNan(SW_PD,3) - meanNoNan(SW_PI,3); 
%SW_estimated_forcing_total2 = meanNoNan(SW_PD - SW_PI,3); 
dat_modis = meanNoNan(SW_estimated_forcing_total,3);
%dat_modis = meanNoNan(SW_estimated_forcing_total2,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-20 10]);

  if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_indirect_estimate_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_indirect_estimate_region_mean = meanNoNan(SW_indirect_estimate_region(:),1);


var_UM = 'Actual indirect change (forcing) in surface SW';
dat_modis = meanNoNan(indirect_ALL,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-20 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 



%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_indirect_actual_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_indirect_actual_region_mean = meanNoNan(SW_indirect_actual_region(:),1)

SW_estimate_regional_mean_bias = 100* ( 1 - SW_indirect_estimate_region_mean./SW_indirect_actual_region_mean);


var_UM = 'Ratio estimated vs actual change';
dat_modis = meanNoNan(SW_estimated_forcing_total,3) ./ meanNoNan(indirect_ALL,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 2]);

var_UM = 'Percentage difference (estimated vs actual change)';
SW_estimated_forcing = SW_PD-SW_PI;
dat_modis = 100*(meanNoNan(SW_estimated_forcing,3) ./ meanNoNan(indirect_ALL,3) - 1);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning2==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

var_UM = 'Change in surface SW due to CF';
SW_estimated_forcing_cf = SW_cf - SW_PI;
SW_estimated_forcing_cf_timemean = meanNoNan(SW_estimated_forcing_cf,3);
dat_modis = SW_estimated_forcing_cf_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 
    
    
 


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_cf_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_cf_region_mean = meanNoNan(SW_cf_region(:),1);


var_UM = 'Change in surface SW due to LWP';
SW_estimated_forcing_lwp = SW_lwp - SW_PI;
SW_estimated_forcing_lwp_timemean = meanNoNan(SW_estimated_forcing_lwp,3);
dat_modis = SW_estimated_forcing_lwp_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_lwp_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_lwp_region_mean = meanNoNan(SW_lwp_region(:),1);


var_UM = 'Change in surface SW due to Nd';
SW_estimated_forcing_Nd = SW_Nd - SW_PI;
SW_estimated_forcing_Nd_timemean = meanNoNan(SW_estimated_forcing_Nd,3);
dat_modis = SW_estimated_forcing_Nd_timemean; 
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 

%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_Nd_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_Nd_region_mean = meanNoNan(SW_Nd_region(:),1);

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
caxis([-20 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_indirect_estimate_linear_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_indirect_estimate_linear_region_mean = meanNoNan(SW_indirect_estimate_linear_region(:),1)

SW_estimate_linear_regional_mean_bias = 100* ( 1 - SW_indirect_estimate_linear_region_mean./SW_indirect_actual_region_mean);



var_UM = 'Percentage bias estimated linear sum vs actual change';
dat_modis =  100*(1 - meanNoNan(SW_estimated_forcing_total_linear,3) ./ meanNoNan(indirect_ALL,3));
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning2==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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

if isave_plots_partioning2==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

%regional mean
SW_estimate_regional_linear_mean_bias_vs_total_estimated = 100* ( 1 - SW_indirect_estimate_linear_region_mean./SW_indirect_estimate_region_mean);



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
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 




var_UM = 'Percentage of abs forcing from \DeltaLWP using PI baseline';
dat_modis =  100*abs(SW_estimated_forcing_lwp_timemean) ./ SW_tot_linear_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 




var_UM = 'Percentage of abs forcing from \DeltaNd using PI baseline';
dat_modis =  100*abs(SW_estimated_forcing_Nd_timemean) ./ SW_tot_linear_timemean_abs;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 

    
return %end the script here in case just want the plots for the figure.



%% ---- Contributions from CF, LWP and Nd using PD with just one for PI switched.
iplot_PI_baseline=0;
if iplot_PI_baseline==1
    
if close_all_figs==1
    close all
    clear gca    
end

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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_cf_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_cf_region_mean_PD = meanNoNan(SW_cf_region(:),1);


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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_lwp_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_lwp_region_mean_PD = meanNoNan(SW_lwp_region_PD(:),1);


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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_Nd_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_Nd_region_mean_PD = meanNoNan(SW_Nd_region_PD(:),1);


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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_indirect_estimate_linear_region_PD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_indirect_estimate_linear_region_mean_PD = meanNoNan(SW_indirect_estimate_linear_region_PD(:),1)

SW_estimate_linear_regional_mean_bias_PD = 100* ( 1 - SW_indirect_estimate_linear_region_mean_PD./SW_indirect_actual_region_mean);


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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

%regional mean
SW_estimate_regional_linear_mean_bias_vs_total_estimated_PD = 100* ( 1 - SW_indirect_estimate_linear_region_mean_PD./SW_indirect_estimate_region_mean);

end

%% Do plots of the paritioning results for PI+PD average - not needed for 2d histos
iplot_PI_baseline=0;
if iplot_PI_baseline==1
    
isave_plots_partioning=0;  %used in draft paper
isave_plots_partioning2=0; %not used in draft at the moment
clear gca

isave_plot=0;


iplot_plots_that_will_be_same=1;

if iplot_plots_that_will_be_same==1
    
    var_UM = 'Total estimated change in surface SW';
    SW_estimated_forcing_total = meanNoNan(SW_PD,3) - meanNoNan(SW_PI,3);
    dat_modis = meanNoNan(SW_estimated_forcing_total,3);
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-20 10]);
    
    if isave_plots_partioning==1
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %        close(gcf);
    end
    
    
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_indirect_estimate_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_indirect_estimate_region_mean = meanNoNan(SW_indirect_estimate_region(:),1);
    
    
    var_UM = 'Actual indirect change (forcing) in surface SW';
    dat_modis = meanNoNan(indirect_ALL,3);
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-20 10]);
    
    if isave_plots_partioning==1
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %        close(gcf);
    end
    
    
    
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
    [iregion_lin,iregion_lin_edges,SW_indirect_actual_region]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis);
    SW_indirect_actual_region_mean = meanNoNan(SW_indirect_actual_region(:),1)
    
    SW_estimate_regional_mean_bias = 100* ( 1 - SW_indirect_estimate_region_mean./SW_indirect_actual_region_mean);
    
    
    var_UM = 'Ratio estimated vs actual change';
    dat_modis = meanNoNan(SW_estimated_forcing_total,3) ./ meanNoNan(indirect_ALL,3);
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 2]);
    
    var_UM = 'Percentage difference (estimated vs actual change)';
    SW_estimated_forcing = SW_PD-SW_PI;
    dat_modis = 100*(meanNoNan(SW_estimated_forcing,3) ./ meanNoNan(indirect_ALL,3) - 1);
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]);
    
    if isave_plots_partioning2==1
        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
        clear opts
        %    opts.iplot_png=1;  %Don't set this since this creates a PNG using
        %    Matlab's driver - but better to convert from eps since otherwise font
        %    is too small
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
    end
    
end

var_UM = 'Change in surface SW due to CF PI+PD average';

SW_estimated_forcing_cf_PIPD_timemean = ( meanNoNan(SW_estimated_forcing_cf_PD,3) + meanNoNan(SW_estimated_forcing_cf,3) )./2;
dat_modis = SW_estimated_forcing_cf_PIPD_timemean;

%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_cf_region_PIPD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_cf_PIPD_region_mean = meanNoNan(SW_cf_region_PIPD(:),1);


var_UM = 'Change in surface SW due to LWP';
SW_estimated_forcing_lwp_PIPD_timemean = ( meanNoNan(SW_estimated_forcing_lwp_PD,3) + meanNoNan(SW_estimated_forcing_lwp,3) )./2;
dat_modis = SW_estimated_forcing_lwp_PIPD_timemean;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_lwp_region_PIPD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_lwp_PIPD_region_mean = meanNoNan(SW_lwp_region_PIPD(:),1);


var_UM = 'Change in surface SW due to Nd';
SW_estimated_forcing_Nd_PIPD_timemean = ( meanNoNan(SW_estimated_forcing_Nd_PD,3) + meanNoNan(SW_estimated_forcing_Nd,3) )./2;
dat_modis = SW_estimated_forcing_Nd_PIPD_timemean; 
%dat_modis = SW_estimated_forcing_Nd_timemean; 
%dat_modis = SW_estimated_forcing_Nd_PD_timemean; 
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 

%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_Nd_region_PIPD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_Nd_PIPD_region_mean = meanNoNan(SW_Nd_region_PIPD(:),1);

% ---- IMPORTANT - here it is important to do the meanNoNan averaging first
% since the CF estimate will have a different number of NaNs than the lwp
% and Nd estimates (since for situations with zero CF in the PI, but some
% CF in PD we are using the PD Nd and LWP values).
var_UM = 'Estimated forcing from linear sum';
SW_estimated_forcing_total_linear_PIPD = SW_estimated_forcing_cf_timemean + SW_estimated_forcing_lwp_PIPD_timemean + SW_estimated_forcing_Nd_PIPD_timemean;
dat_modis =  SW_estimated_forcing_total_linear_PIPD;
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-20 10]);

    if isave_plots_partioning==1       
        savename=[savedir_date titlenam_driver];        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
%        close(gcf);
    end 


%restrict to the region of interest
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin,iregion_lin_edges,SW_indirect_estimate_linear_region_PIPD]=get_lat_lon_irregular_with_time(1,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,dat_modis); 
SW_indirect_estimate_linear_region_mean_PIPD = meanNoNan(SW_indirect_estimate_linear_region_PIPD(:),1);

SW_estimate_linear_regional_mean_bias_PIPD = 100* ( 1 - SW_indirect_estimate_linear_region_mean_PIPD./SW_indirect_actual_region_mean);



var_UM = 'Percentage bias estimated linear sum vs actual change';
dat_modis =  100*(1 - meanNoNan(SW_estimated_forcing_total_linear_PIPD,3) ./ meanNoNan(indirect_ALL,3));
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning2==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end


%I.e., linear sum vs total change, both using SW estimates from the
%calculation.
var_UM = 'Percentage bias estimated linear sum vs estimated total change';
dat_modis =  100*(1 - meanNoNan(SW_estimated_forcing_total_linear_PIPD,3) ./ meanNoNan(SW_estimated_forcing,3));

%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

if isave_plots_partioning2==1
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_' titlenam_driver];
    clear opts
%    opts.iplot_png=1;  %Don't set this since this creates a PNG using
%    Matlab's driver - but better to convert from eps since otherwise font
%    is too small
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

%regional mean
SW_estimate_regional_linear_mean_bias_vs_total_estimated_PIPD = 100* ( 1 - SW_indirect_estimate_linear_region_mean_PIPD./SW_indirect_estimate_region_mean);

end

%% Plots for PD baseline showing the proportion of forcing attributable to each variable
iplot_PI_baseline=0;
if iplot_PI_baseline==1
    
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
end

end

%% Do the calculations to pick out the regions of large positive CF change
% between PI and PD
iplot_large_CF=0;
if iplot_large_CF==1
    
if close_all_figs==1
    close all
    clear gca    
end

%dlowcf_timemean is the change in low CF between PI and PD

%Using the same routine as for the lat lon region, but for dcf. Will
%replace the lon values with dummy values and ranges that are always satisfied
dummy_gcm_Plon2D_UM = ones(size(dlowcf_timemean));
dummy_LON_thresh=[0 2];
%dcf_thresh=[0.05 9e9];
dcf_thresh=[0.03 9e9];
[iregion_lin_dcf,iregion_lin_edges_dcf,cf_low_dcf]=get_lat_lon_irregular_with_time(nT_main,dcf_thresh,dummy_LON_thresh,dlowcf_timemean,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,low_CF_PI_ALL(:,:,1:end));
[iregion_lin_dcf,iregion_lin_edges_dcf,cf_low_dcf_PD]=get_lat_lon_irregular_with_time(nT_main,dcf_thresh,dummy_LON_thresh,dlowcf_timemean,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,low_CF_PD_ALL(:,:,1:end));

[iregion_lin_dcf,iregion_lin_edges_dcf,cf_mid_dcf]=get_lat_lon_irregular_with_time(nT_main,dcf_thresh,dummy_LON_thresh,dlowcf_timemean,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,mid_CF_PI_ALL(:,:,1:end));
[iregion_lin_dcf,iregion_lin_edges_dcf,cf_mid_dcf_PD]=get_lat_lon_irregular_with_time(nT_main,dcf_thresh,dummy_LON_thresh,dlowcf_timemean,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,mid_CF_PD_ALL(:,:,1:end));

[iregion_lin_dcf,iregion_lin_edges_dcf,cf_high_dcf]=get_lat_lon_irregular_with_time(nT_main,dcf_thresh,dummy_LON_thresh,dlowcf_timemean,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,high_CF_PI_ALL(:,:,1:end));
[iregion_lin_dcf,iregion_lin_edges_dcf,cf_high_dcf_PD]=get_lat_lon_irregular_with_time(nT_main,dcf_thresh,dummy_LON_thresh,dlowcf_timemean,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,dummy_gcm_Plon2D_UM,high_CF_PD_ALL(:,:,1:end));



%Now restrict to the box region too
LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;
[iregion_lin_dcf_box,iregion_lin_edges_dcf_box,cf_low_dcf_box]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,cf_low_dcf); 
[iregion_lin_dcf_box,iregion_lin_edges_dcf_box,cf_low_dcf_box_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,cf_low_dcf_PD);     

[iregion_lin_dcf_box,iregion_lin_edges_dcf_box,cf_mid_dcf_box]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,cf_mid_dcf); 
[iregion_lin_dcf_box,iregion_lin_edges_dcf_box,cf_mid_dcf_box_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,cf_mid_dcf_PD);     

[iregion_lin_dcf_box,iregion_lin_edges_dcf_box,cf_high_dcf_box]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,cf_high_dcf); 
[iregion_lin_dcf_box,iregion_lin_edges_dcf_box,cf_high_dcf_box_PD]=get_lat_lon_irregular_with_time(nT_main,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,cf_high_dcf_PD);     


var_UM = '';
dat_modis =  meanNoNan(cf_low_dcf_box,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([0 100]);

var_UM = '';
dat_modis =  meanNoNan(cf_low_dcf_box_PD,3);
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis([0 100]);

end

%% 2D mean contribution to forcing from different combinations of PI and PD low cloud fraction  - Regions for high dcf only
% 
iplot_large_CF=0;
if iplot_large_CF==1
    
    clear gca
    isave_mean_2D_PDF=0;
    
    
% Plot the forcing vs the cloud fraction in states 1 and 2 (1=clear -
% included to get zero CFs).

    %inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    %inds_PI_clear_low = unique( [inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    
    %I think that the way it is done below will allow non-clear/low points
    %to be included for one of the states
    %inds_PI_clear_low = unique( [inds_PI{1}; inds_PI{2}; inds_PD{1}; inds_PD{2}] ); %Both states 1 and 2 (clear and low-only)
    
    %Whereas using intersect here should make it so that e.g. PI needs to
    %be in the clear state AND PD in the clear state (but using the 4
    %combinations :- clear + clear, clear + low, etc.)
    inds_CC = intersect( inds_PI{1}, inds_PD{1} );
    inds_CL = intersect( inds_PI{1}, inds_PD{2} );
    inds_LC = intersect( inds_PI{2}, inds_PD{1} );
    inds_LL = intersect( inds_PI{2}, inds_PD{2} );
    
    inds_PI_clear_low = unique( [inds_CC; inds_CL; inds_LC; inds_LL] );
        
    X_driver = cf_low_dcf_box(inds_PI_clear_low);
    Y_driver = cf_low_dcf_box_PD(inds_PI_clear_low);    
%    Z_driver = dlowcf_ALL(inds_PI_clear_low); %forcing(inds_PI_clear_low);
    Z_driver = forcing(inds_PI_clear_low);    
    xlabelstr='Low cloud fraction Pre-Industrial';
    ylabelstr = 'Low cloud fraction Present Day';
    
%plot the mean forcing in each cloud fraction bin
    %run 2D histo template script :-
    DRIVER_template_2D_mean_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-100 100]);
    title('Mean forcing (W m^{-2})');
    %shading faceted; %this adds grid lines to the plot - but gives NaNs a
    %colour...
    label_text_pcolor; %Script to add text labels of the numbers for each block
    save_str = 'forcing_mean_2D_PDF_cf_PI_PD';
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  
        
%plot the overaell forcing (mean * frequency) in each cloud fraction bin    - i.e., these sum to the total forcing 
    %N.B. - the difference here is that the following scripts sets this
    %flag :- ioverall_contribution_mean_2D = 1;
    DRIVER_template_2D_overall_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis([-5 5]);    
    %caxis([-2.5 2.5]);   
    %shading faceted; %this adds grid lines to the plot
    label_text_pcolor; %Script to add text labels of the numbers for each block
    title('Overall forcing (W m^{-2})');
    save_str = 'forcing_overall_2D_PDF_dcf_region'; 
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] ); 
    
end

%% 2D plots of forcing vs cloud state PI and PD for regions of high CF
iplot_large_CF=0;
if iplot_large_CF==1
    
    X_driver = cloud_states_PI(iregion_lin_dcf_box);
    Y_driver = cloud_states_PD(iregion_lin_dcf_box);    
    Z_driver = forcing(iregion_lin_dcf_box);
    xlabelstr='Cloud State Pre-Industrial';
    ylabelstr = 'Cloud State Present Day';    
    
    DRIVER_template_2D_overall_SW_forcing_vs_cloudstate_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis([-2 2]);    
    %caxis([-2.5 2.5]);       
    title('Overall forcing (W m^{-2})');
    label_text_pcolor; %Script to add text labels of the numbers for each block
    set(gca,'xlim',[0.5 8.5]);
    set(gca,'ylim',[0.5 8.5]);    
    save_str = 'forcing_overall_2D_PDF_vs_cloud_states_PI_PD';
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  
    
        
    
end

%% Do the same, but for rain rates
iplot_large_CF=0;
if iplot_large_CF==1
    
if close_all_figs==1
    close all
    clear gca    
end

%N.B. - don't need to worry about restricting the rain rates themselves to
%the box region, etc. since the X_driver and Y_driver data have already
%been restricted.

    X_driver = cf_low_dcf_box(inds_PI_clear_low);
    Y_driver = cf_low_dcf_box_PD(inds_PI_clear_low);    
%    Z_driver = dlowcf_ALL(inds_PI_clear_low); %forcing(inds_PI_clear_low);
    Z_driver = 24*3600*(surf_rain_rate_PD_ALL(inds_PI_clear_low) - surf_rain_rate_PI_ALL(inds_PI_clear_low));    
    xlabelstr='Low cloud fraction Pre-Industrial';
    ylabelstr = 'Low cloud fraction Present Day';
    
%plot the mean forcing in each cloud fraction bin
    %run 2D histo template script :-
    DRIVER_template_2D_mean_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(24*[-0.07 0.07]);
    caxis(24*[-0.02 0.02]);
    title('Bin mean surface rain rate change (mm day^{-1})');        
    label_text_pcolor; %Script to add text labels of the numbers for each block
    save_str = 'rain_change_mean_2D_PDF_cf_PI_PD';
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  
        
%plot the overaell forcing (mean * frequency) in each cloud fraction bin    - i.e., these sum to the total forcing 
    %N.B. - the difference here is that the following scripts sets this
    %flag :- ioverall_contribution_mean_2D = 1;
    DRIVER_template_2D_overall_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis(24*[-4 4]*1e-4);    
    %caxis([-2.5 2.5]);       
    title('Bin contribution to surface rain rate change (mm day^{-1})');
    label_text_pcolor; %Script to add text labels of the numbers for each block    
    save_str = 'rain_rate_overall_2D_PDF_dcf'; 
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  

    inan=find(isnan(X_driver)==1);
    surf_rain_rate_dcf_box=surf_rain_rate_PI_ALL(inds_PI_clear_low); surf_rain_rate_dcf_box(inan)=NaN;
    surf_rain_rate_dcf_box_PD=surf_rain_rate_PD_ALL(inds_PI_clear_low); surf_rain_rate_dcf_box_PD(inan)=NaN;
    
    
end

%% Rain rates divided by cloud fraciton (in-cloud rain rates)
%N.B. - don't need to worry about restricting the rain rates themselves to
%the box region, etc. since the X_driver and Y_driver data have already
%been restricted.
iplot_large_CF=0;
if iplot_large_CF==1
    
in_cloud_rain_box = surf_rain_rate_PI_ALL./cf_low_dcf_box; in_cloud_rain_box(cf_low_dcf_box<cf_thresh)=NaN;
in_cloud_rain_box_PD = surf_rain_rate_PD_ALL./cf_low_dcf_box_PD; in_cloud_rain_box_PD(cf_low_dcf_box_PD<cf_thresh)=NaN;

    X_driver = cf_low_dcf_box(inds_PI_clear_low);
    Y_driver = cf_low_dcf_box_PD(inds_PI_clear_low);    
%    Z_driver = dlowcf_ALL(inds_PI_clear_low); %forcing(inds_PI_clear_low);
    Z_driver = 24*3600*(in_cloud_rain_box_PD(inds_PI_clear_low) - in_cloud_rain_box(inds_PI_clear_low));    
    xlabelstr='Low cloud fraction Pre-Industrial';
    ylabelstr = 'Low cloud fraction Present Day';
    
%plot the mean forcing in each cloud fraction bin
    %run 2D histo template script :-
    DRIVER_template_2D_mean_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(24*[-1.4 1.4]);
    title('Bin mean in-cloud surface rain rate change (mm day^{-1})');
    label_text_pcolor; %Script to add text labels of the numbers for each block    
    save_str = 'forcing_mean_2D_PDF_cf_PI_PD';
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  
        
%plot the overaell forcing (mean * frequency) in each cloud fraction bin    - i.e., these sum to the total forcing 
    %N.B. - the difference here is that the following scripts sets this
    %flag :- ioverall_contribution_mean_2D = 1;
    DRIVER_template_2D_overall_SW_forcing_vs_CF_UM   
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    %caxis([-3.1 3.1]);
    caxis(24*[-0.02 0.02]);    
    %caxis([-2.5 2.5]);       
    title('Bin contribution to in-cloud surface rain rate change (mm day^{-1})');
    label_text_pcolor; %Script to add text labels of the numbers for each block
    save_str = 'rain_rate_overall_2D_PDF_dcf'; 
    savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
    UM_save_plot(gcf,isave_mean_2D_PDF,savedir,[save_str] );  

    inan=find(isnan(X_driver)==1);
    in_cloud_surf_rain_rate_dcf_box=in_cloud_rain_box(inds_PI_clear_low); in_cloud_surf_rain_rate_dcf_box(inan)=NaN;
    in_cloud_surf_rain_rate_dcf_box_PD=in_cloud_rain_box_PD(inds_PI_clear_low); in_cloud_surf_rain_rate_dcf_box_PD(inan)=NaN;
    
end    


irun_extra_code=0;
if irun_extra_code==1


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

%% to plot the region box can run this code

ioverride_box_colour=1;

LAT_val = LAT_val_DRIVER2;
LON_val = LON_val_DRIVER2;

itext_in_box=0;
box_lwidth=4; %linewidth of box

irotated_pole_box=0;

imap=1;

plot_box_on_map
        
%% plots that probably don't need any more
iplot_extras=0;
if iplot_extras==1
    

    
%% bar charts of Contribution from cloud states
    


    
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
%         saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        %    close(gcf);
    end
    
    
     
    
    % ------------------
    
%% Forcing and 2D histograms by cloud fraction
clear gca
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
%     saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
%     saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
%     saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
%     saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
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
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%    close(gcf);
end 

    

    
end


end