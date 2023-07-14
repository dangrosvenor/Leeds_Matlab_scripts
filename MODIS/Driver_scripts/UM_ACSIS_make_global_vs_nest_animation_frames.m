%Add variables to climits are set here too :-
%   UM_var_defs.m

%runs this plotting script - lat/lon for map also set here :-
  %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
%FOR PDFs, etc. it is set below
  
load_type = 'mat';
load_type = 'merged netCDF';

cloud_input = 'UM'; %Use the usual UM low, mid and high cloud fractions.
%cloud_input = 'CALIPSO'; %Use the UM COSP CALIPSO values
%cloud_input = 'MODIS'; %Use the COSP MODIS values

%[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

amsre_data = 1;

ical_data=1; %Whether are using CALIPSO data - load using script :- read_calipso_monthly_night_IPSL.m for average values
    %first. This puts the data into these fields :- 
    % cllcalipso_monthly_AVERAGE, clmcalipso_monthly_AVERAGE, clhcalipso_monthly_AVERAGE
    % Is average of day and night values for now (have separate ones too).

iceres_data=1; %whether to load CERES data
    
%um_case_PD = 'u-av503';
um_case_nest = 'u-at459';
um_case_global = 'u-at459_glm';

time_format_str=' UTC';

LAT_val_DRIVER2 = [30 45]; LON_val_DRIVER2 = [-60 -10]; %N Atlantic region of high negative forcing ocean only
LAT_val_DRIVER2 = [30 50]; LON_val_DRIVER2 = [-60 -10]; %N Atlantic region of high negative forcing ocean only

%  /home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-au652/SW_down_surf//umglaa_SW_down_surf_native res_ALL.mat    %PI
%  /home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-au536/SW_down_surf//umglaa_SW_down_surf_native res_ALL.mat     %PD

iplot_global=1; %map of SW forcing

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


var_UM_DRIVER = 'LWP';
%var_UM_DRIVER = 'accum_number_ukca';
%var_UM_DRIVER = 'SW_down_surf';
%var_UM_DRIVER = 'Nd_lwc_weighted_UKCA';


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
clear time_choice; time_choice.time_range = [datenum('01-Jan-0000') datenum('31-Dec-9999')];

%Specific days :-
%clear time_choice; time_choice.time_specific = datenum('06-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('07-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('10-Aug-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Jul-2009'); time_choice.find_nearest=1;
%clear time_choice; time_choice.time_specific = datenum('31-Dec-2009'); time_choice.find_nearest=1;

%% maps of LWP for nest and global

scrsz=get(0,'ScreenSize');
posit=[9 60 scrsz(3) scrsz(4)];
%figure('position',posit);
figure('position',scrsz);




if iplot_global==1
    %subplot(1,3,1);
    var_UM = var_UM_DRIVER;

% PI run    
%    um_case='u-au652'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    um_case=um_case_global; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions    
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];

    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

    

array_in=[]; %just test to get the indices for now.
dim=NaN; %don't need dim if just getting the indices
[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);
time_round = datestr(dat_global.time_ALL(time_inds(1)));


    
    gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
    gcm_Plat2D_edges_UM = dat_global.gcm_Plat2D_edges_UM;

    gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
    %    if maxALL(gcm_Plon2D_UM) < 181
    if irestrict_domain_DRIVER==0
        i180 = find(gcm_Plon2D_UM<0);
        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180)+360;
    end
    [gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);

    
    [LWP_global_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
    
    %restrict to the region of interest
    LAT_val = LAT_val_DRIVER2;
    LON_val = LON_val_DRIVER2;
       
    [iregion_lin,iregion_lin_edges,LWP_global_ALL_region]=get_lat_lon_irregular_with_time(nT,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,LWP_global_ALL); 

    
    
    var_UM = 'LWP'; 
 
    um_case=um_case_nest; pole_lat=45.0; pole_lon=145.0; run_type = 'nested'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];    
    dat_nest = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon);  
    [LWP_nest_ALL,nT] = UM_get_time_data_mat_nc(dat_nest,time_inds,load_type,gcm_Plat2D_UM);
   
    
    for it_anim=1:size(LWP_nest_ALL,3);
        UM_ACSIS_global_nest_animation_SUBPLOT_commands
    end
    
    
    



end

