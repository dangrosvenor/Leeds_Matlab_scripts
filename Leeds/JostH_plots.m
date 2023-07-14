%SO timseries for Jost Heizenberg work

%Data in :-
%/home/disk/eos15/d.grosvenor/UM/UKESM/UKESM_ACSIS_nudged_run_James_Keeble/u-bz785/

%Based on UM_Iceland_plot_various_based_on_ACSIS_script.m


%Add variables to climits are set here too :-
%   UM_var_defs.m

%runs this plotting script - lat/lon for map also set here. FOR PDFs, etc. it is set below :-
  %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
  

UM_ACSIS_SW_vs_cloud_properties_global_DEFAULTS

savedir_date=['/home/disk/eos1/d.grosvenor/UM/JostH/plots_' datestr(now,30) '/'];
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
run_set = 'Blending option=3, volcano starting 12 UTC';

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
        um_case_PI = 'u-ch765'; run_type_DRIVER = 'nested ignore lat lon'; %volcano OFF
        um_case_PD = 'u-ch764'; run_type_DRIVER = 'nested ignore lat lon'; %volcano ON   
        Nd_var_str='UKCA';  
        Nd_var_str2='_total_column_to_zdomain_top';  
        
        
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
var_UM_DRIVER = 'N_Dgt3'; %use the file, but only for loading the time

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
clear time_choice; time_choice.time_range = [datenum('20-Dec-2020') datenum('24-Dec-2020 13:00')]; %
clear time_choice; time_choice.time_range = [datenum('20-Dec-2020') datenum('22-Dec-2020 12:00')]; %
clear time_choice; time_choice.time_range = [datenum('01-Jan-0000') datenum('31-Dec-3000 12:00')]; %

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
            um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = 'nested'; icoarse=0; ivar_dir=1; %nested UKCA
            dirUM = ['/home/disk/eos15/d.grosvenor/UM/Hawaii/' um_case];
            
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

%save_time_lat_file = [dirUM '/save_time_lat.mat'];
%save(save_time_lat_file,'-V7.3');
    


%%
      load_type = 'merged netCDF';
        var_UM = 'Land_mask'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
        um_case='u-bf666'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global_mask = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
        
        

%% N D>3nm


time_inds = [1:12];

run_type_DRIVER = 'global';
load_type = 'merged netCDF';

% Load in the N for D>3nm
    var_UM = 'N_Dgt3';
    %var_UM = 'SO2_perm3_total_column_to_zdomain_top';
    %var_UM2 = 'SO2_perm3';
    um_case='u-bz785'; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/UKESM/UKESM_ACSIS_nudged_run_James_Keeble/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon); 
    
    
%Need the average over months 10 to 4 for     
iz=1; %Choose the lowest model level for now.
[dat_surf_time_mean,N,std_dev] = meanNoNan(dat_global.dat([1:4 10:12],iz,:,:),1);
dat_surf_time_mean(dat_global_mask.dat==1) = NaN;

%Now need to plot zonal means for lat bands : -75 to -65, etc. to -5
lats = [-75:10:-5];
%dat = dat_global.dat([1:4 10:12],iz,:,:);

clear dat_means dat_N dat_std xline yline
xline(1) = lats(1); yline(1)=0;
for i=1:length(lats)-1
   imean = find( dat_global.gcm_Plat2D_UM > lats(i) & dat_global.gcm_Plat2D_UM <= lats(i+1) );   
   [dat_means(i),dat_N(i),dat_std(i)] = meanNoNan(dat_surf_time_mean(imean),1);    
   xline(end+1) = lats(i); yline(end+1) = dat_means(i);
   xline(end+1) = lats(i+1); yline(end+1) = dat_means(i);
   xline(end+1) = lats(i+1); yline(end+1) = 0;
end
%xline(end+1) = lats(end); yline(end+1) = 0;

figure
set(gcf,'color','w');
mid_points = 0.5*(lats(1:end-1) + lats(2:end));
%plot(mid_points,dat_means/1e6,'bo');
%set(gca,'xlim',[-75 0]);
line(xline,yline/1e6,'linewidth',3);
ylabel('Mean N_{D>3nm}');
xlabel('Latitude');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-75 0]);

hold on
herr=errorbarYY('vert',mid_points,dat_means/1e6,dat_std/1e6,gca,'b','o',2,0.01);
%errorbar(mid_points,dat_means/1e6,dat_std/1e6);
set(gca,'ylim',[0 800]);   
    
    
    
%% N D>10


time_inds = [1:12];

run_type_DRIVER = 'global';
load_type = 'merged netCDF';

% Load in the N for D>3nm
    var_UM = 'N_Dgt10';
    %var_UM = 'SO2_perm3_total_column_to_zdomain_top';
    %var_UM2 = 'SO2_perm3';
    um_case='u-bz785'; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %
    dirUM = ['/home/disk/eos15/d.grosvenor/UM/UKESM/UKESM_ACSIS_nudged_run_James_Keeble/' um_case];    
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon); 
    
    
%Need the average over months 10 to 4 for     
iz=1; %Choose the lowest model level for now.
[dat_surf_time_mean,N,std_dev] = meanNoNan(dat_global.dat([1:4 10:12],iz,:,:),1);
dat_surf_time_mean(dat_global_mask.dat==1) = NaN;

%Now need to plot zonal means for lat bands : -75 to -65, etc. to -5
lats = [-75:10:-5];
%dat = dat_global.dat([1:4 10:12],iz,:,:);

clear dat_means dat_N dat_std xline yline
xline(1) = lats(1); yline(1)=0;
for i=1:length(lats)-1
   imean = find( dat_global.gcm_Plat2D_UM > lats(i) & dat_global.gcm_Plat2D_UM <= lats(i+1) );   
   [dat_means(i),dat_N(i),dat_std(i)] = meanNoNan(dat_surf_time_mean(imean),1);    
   xline(end+1) = lats(i); yline(end+1) = dat_means(i);
   xline(end+1) = lats(i+1); yline(end+1) = dat_means(i);
   xline(end+1) = lats(i+1); yline(end+1) = 0;
end
%xline(end+1) = lats(end); yline(end+1) = 0;

figure
set(gcf,'color','w');
mid_points = 0.5*(lats(1:end-1) + lats(2:end));
%plot(mid_points,dat_means/1e6,'bo');
%set(gca,'xlim',[-75 0]);
line(xline,yline/1e6,'linewidth',3);
ylabel('Mean N_{D>10nm} (cm^{-3})');
xlabel('Latitude');
fontsize_figure(gcf,gca,18);
set(gca,'xlim',[-75 0]);

hold on
herr=errorbarYY('vert',mid_points,dat_means/1e6,dat_std/1e6,gca,'b','o',2,0.01);
%errorbar(mid_points,dat_means/1e6,dat_std/1e6);
set(gca,'ylim',[0 800]);  