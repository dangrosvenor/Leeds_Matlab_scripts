% Timeseries plot for UM

%Choose runs here :- UM_case_select_runs

%Timeseries loaded (using load_UM_timser) and added together here :- UM_add_timeseries

isave_plot_driver=0;

i_restricted_domain=1; %e.g. to add _restricted_domain_.pp.nc in etc. - have added this directly in UM_add_timeseries

inormalise_initial=0; %Normalise by dividing by first value (relative diff)
inormalise_initial_absolute=1; %Normalise by subtracting first value
iadd_umid_label=1; %Add the Run ID (e.g. u-ar365) to the legend labels

iplot_flux_contribution=0; %Plot the domain mean contribution from the flux rather than the actual domain mean
iplot_flux_increments=0; %Calculate the increments from the flux values and compare them to the actual change (for conserved quantities)

clear ztop_int

%var_tim = 'droplet_number_total_column_to_z1500';
% var_tim = 'Accum mode aerosol number column integrated z1500m budget';
% var_tim = 'coarse_number_total_column_to_z1500';
%var_tim = 'coarse_mass_total_column_to_z1500';
%var_tim = 'accum_mass_total_column_to_z1500';
%var_tim = 'act_mass_liq_total_column_to_z1500';
%var_tim = 'act_mass_rain_total_column_to_z1500';

%var_tim = 'Coarse mode aerosol number column integrated z1500m budget';

% var_tim = 'Total aerosol number column integrated z1500m budget'; %N.B.-this includes droplet number
% var_tim = 'Total aerosol mass column integrated z1500m budget'; %Aerosol mass inc. activated

%var_tim = 'Total aerosol mass no act column integrated z1500m budget'; %Aerosol mass without activated aerosol, i.e. intestitial only

% --- z1500 to domain top ---
%var_tim = 'Total aerosol number column integrated budget z1500m to domain top';

% var_tim = 'droplet_number_column_integrated';
% var_tim = 'Total aerosol number column integrated budget';  %For other ones where are adding stuff together, etc. need an extra case later on

% --- to z3000 ---

var_tim = 'Total aerosol number column integrated z3000m budget';
var_tim = 'accum_number_total_column_to_z3000';  %If it is just a single standard variable then jsut have to set its name and the title below
%var_tim = 'droplet_number_total_column_to_z3000'; %but the normalization won't work
var_tim = 'coarse_number_total_column_to_z3000'; %but the normalization won't work
%var_tim = 'Droplet number column integrated to specified z'; ztop_int=3000; %normalisation will work with this one.

%var_tim = 'Accumulation mode aerosol mass column integrated to specified z'; ztop_int=3000; %normalisation will work with this one.
%var_tim = 'Aitken mode aerosol mass column integrated to specified z'; ztop_int=3000; %normalisation will work with this one.
%var_tim = 'Aitken mode aerosol number column integrated to specified z'; ztop_int=3000; %normalisation will work with this one.
%var_tim = 'Coarse mode aerosol mass column integrated to specified z'; ztop_int=3000; %normalisation will work with this one.
%var_tim = 'Activated aerosol mass column integrated to specified z'; ztop_int=3000; %normalisation will work with this one.
%var_tim = 'Air mass column integrated to specified z'; ztop_int=3000;
var_tim = 'Total aerosol mass column integrated to specified z'; ztop_int=3000;


%var_tim = 'qL column integrated to specified z'; ztop_int=3000; %normalisation will work with this one.

% --- to domain top ---
%var_tim = 'Air mass column integrated to specified z'; ztop_int=1e9;

%var_tim = 'Total aerosol number column integrated to domain top budget';
%var_tim = 'Coarse mode aerosol number column integrated to domain top budget'; %
%var_tim = 'Accum mode aerosol number column integrated to domain top budget'; %

%var_tim = 'Droplet number column integrated to domain top budget'; %

%var_tim = 'Total aerosol mass column integrated to domain top budget'; %Aerosol mass inc. activated
%var_tim = 'Coarse mode aerosol mass column integrated to domain top budget'; %
%var_tim = 'Accum mode aerosol mass column integrated to domain top
%budget'; %

%var_tim = 'Activated in-cloud aerosol mass column integrated to domain top budget';




savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

iplot_goes=0; %Whether to plot the GOES timeseries
iplot_RHB=0;  %Whether to plot the RHB ship timeseries
iplot_amsre=0;
ismooth_RHB=0;
isave_plot=0;

day_or_night = 'all'; %Used for the calculation of Nd_overall_mean and LWP_overall_mean

% --- Provide the name of the set of runs to process (see :-
%       --- UM_case_select_runs ---
% function)
UM_cases = '12th Nov case, as of May 2016';
UM_cases = '12th Nov case, as of May 2016 adhoc';
UM_cases = '12th Nov case, as of May 2016 adhoc eos10';
UM_cases = '12th Nov case, as of May 2016 adhoc multi-dirUM';
%UM_cases = '12th Nov case, as of May 2016 processing runs multi-dirUM';
%UM_cases = '12th Nov case, as of May 2016 rain OFF multi-dirUM';
UM_cases = '12th Nov case, as of May 2016 processing runs PLOTS multi-dirUM';
UM_cases = '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM';
UM_cases = '12th Nov case, as of May 2017 processing runs surface fluxes, delaying processing, etc.';
%UM_cases = '12th Nov case, cloud scheme test v10.8';

%UM_cases = 'Iceland_9day_runs_Nov2016';

% --- Runs this script to get the filenames :- UM_case_select_RUN  

%% Lat, lon, etc. - used for timeseries?

LAT_val_DRIVER = [-23.5 -16.44]; LON_val_DRIVER = [-85.93 -78.08]; %Smaller region to but out model edges
LAT_val_DRIVER = [-24.5 -15.44]; LON_val_DRIVER = [-86.93 -77.08]; %GOES region for UM comparison xkqk 26thOct POC
LAT_val_DRIVER = [-21 -16.44]; LON_val_DRIVER = [-82 -78.08]; %Much smaller region in NW corner
LAT_val_DRIVER = [-24.5 -15.44]; LON_val_DRIVER = [-84 -77.08]; %Same region as for the AGU PDFs - Trying to match AMSRE and GOES domains
LAT_val_DRIVER = [-21 -15.44]; LON_val_DRIVER = [-86.93 -82]; %Revised PDFs post-AGU. Revised smaller region in NW corner  

% 12th Nov case
%Need to come up with a new corner for the 12th Nov case since this is
%centred at 76W, 20S (not at 82W).
% Domain runs from -78.93 to -73.08 W and -22.70 to -17.28 (based on the edges)

%FULL UM domain for 12th Nov :-
LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER = [-78.93 -73.08];

thresh_LWP_DRIVER = 5; %threshold for UM Nd - sets values with lower LWP to NaN, so that they aren't included in the mean

% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];

%The location of the ship = 20S, 75W
%One degree region around it
%LAT_val_DRIVER = [-20.5 -19.5]; LON_val_DRIVER = [-75.5 -74.5];


%LAT_val_DRIVER = [-20.5 -20.4]; LON_val_DRIVER = [-75.1 -75.0];

% See --- read_GOES_vocals_netcdf_files_MULTI.m ---
%           for reading in muliple files and saving the data in a .mat file
%           26-27th Oct
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20150509T124841.mat';
%  12-14th Nov
%   LAT_val = [-22.70 -17.28]; LON_val = [-78.93 -73.08];
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20150730T092051.mat';

%Using partial domain to avoid boundary issues
%   LAT_val= [-20.5 -17.5]; LON_val = [-78.75 -73.25];
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20150731T024411.mat';
load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151112T064655.mat';
  %Special test case for 1 deg region around the ship to see if GOES
  %compares well with the aircraft :-
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151112T074537.mat';  

thresh_SZA = 64;

%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150727T071106.mat';
%Using full domain
%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150729T083601.mat';
%Using partial domain to avoid boundary issues
%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150729T101619.mat';


%runt his script to get cdan(1:9) and markers(1:3) and pdan(1:3)
LInestyles_etc

%--- run the file to set up the defaults
watervap_defaults
for idat=1:99
    ismooth_x_import(idat)=0;
    ismooth_y_import(idat)=0;
end
idat_driver=1;

clear fileUM xdat_import ydat_import flag 

%--- set some options for this particular plot
graph=0; %graph choice in watervap

xlab='Time (UTC)';
xlab='Time (Local Solar Time)';


if iplot_flux_increments==1 %if running this then are plotting two variables for each run
    idat_driver_inc=2;
else
    idat_driver_inc=1;
end


if exist('ztop_int')
    if ztop_int>999e3
        hgt_str = 'domain_top';
        hgt_str2 = 'domain top';
    else
        hgt_str = num2str(ztop_int);
        hgt_str2 = hgt_str;
    end
else
    hgt_str='';
    hgt_str2='';
end



switch var_tim
     case 'coarse_number_total_column_to_z3000'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated coarse mode number concentration, z<3km (m^{-2})';
    case 'accum_number_total_column_to_z3000'
        ylab='Domain-mean (m^{-2})';
        titlenam = 'Timeseries of domain-mean column integrated accum. mode aerosol number concentration, z<3km';
    case 'droplet_number_total_column_to_z3000'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated droplet number concentration, z<3km (m^{-2})';
    case 'droplet_number_total_column_to_z1500'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated droplet number concentration, z<1.5km (m^{-2})';
    case 'coarse_number_total_column_to_z1500'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated coarse mode number concentration, z<1.5km (m^{-2})';
    case 'coarse_mass_total_column_to_z1500'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated coarse mode mass MR, z<1.5km (kg m^{-2})';
    case 'accum_mass_total_column_to_z1500'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated accumulation mode mass MR, z<1.5km (kg m^{-2})';         
    case 'act_mass_liq_total_column_to_z1500'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated mass MR in droplets, z<1.5km (kg m^{-2})';
    case 'act_mass_rain_total_column_to_z1500'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated mass MR in rain, z<1.5km (m^{-2})';
    case 'Total aerosol number column integrated budget'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated aerosol number concentration';   
    case 'Total aerosol number column integrated z3000m budget'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated aerosol number concentration to z=3000m';   
    case 'Total aerosol number column integrated z1500m budget'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated aerosol number concentration to z=1500m';     

    case 'Total aerosol number column integrated to domain top budget'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated aerosol number concentration to domain top';             
        
    case 'Total aerosol mass column integrated to specified z'
        ylab='Domain-mean (kg m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated total aerosol mass mixing ratio to z=' hgt_str2];
        
     case 'Air mass column integrated to specified z'
        ylab='Domain-mean (kg m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated air mass mixing ratio to z=' hgt_str2];
   
     case 'Aitken mode aerosol number column integrated to specified z'
        ylab='Domain-mean (m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated Aitken mode aerosol num mixing ratio to z=' hgt_str2];
        
     case 'Aitken mode aerosol mass column integrated to specified z'
        ylab='Domain-mean (kg m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated Aitken mode aerosol mass mixing ratio to z=' hgt_str2];
           
    case 'Accumulation mode aerosol mass column integrated to specified z'
        ylab='Domain-mean (kg m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated accum mode aerosol mass mixing ratio to z=' hgt_str2];

    case 'Coarse mode aerosol mass column integrated to specified z'
        ylab='Domain-mean (kg m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated coarse mode aerosol mass mixing ratio to z=' hgt_str2];
        

    case 'Activated aerosol mass column integrated to specified z'
        ylab='Domain-mean (kg m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated act aerosol mass mixing ratio to z=' hgt_str2];
        
        
     case 'Total aerosol mass column integrated z1500m budget'
        ylab='Domain-mean (kg m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated total aerosol mass mixing ratio to z=1500m'; 
     case 'Total aerosol mass column integrated to domain top budget'    
         ylab='Domain-mean (kg m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated aerosol mass mixing ratio to domain top';         
     case 'Coarse mode aerosol mass column integrated to domain top budget'    
         ylab='Domain-mean (kg m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated coarse mode aerosol mass mixing ratio to domain top';  
      case 'Coarse mode aerosol number column integrated to domain top budget'    
         ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated coarse mode aerosol number concentration to domain top';   
      case 'Accum mode aerosol number column integrated to domain top budget'    
         ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated accumulation mode aerosol number concentration to domain top';               
      case 'Droplet number column integrated to domain top budget'    
         ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated droplet number concentration to domain top';                       
    case 'Droplet number column integrated to specified z'
        ylab='Domain-mean (m^{-2})';
        titlenam=['Timeseries of domain-mean of column integrated N_d to z=' hgt_str2];
    case 'qL column integrated to specified z'
        ylab='Domain-mean (kg m kg^{-1})';
        titlenam=['Timeseries of domain-mean of column integrated qL to z=' hgt_str2];        
     case 'Accum mode aerosol mass column integrated to domain top budget'    
         ylab='Domain-mean (kg m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated accumulation mode aerosol mass mixing ratio to domain top';                 
     case 'Activated in-cloud aerosol mass column integrated to domain top budget'
         ylab='Domain-mean (kg m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated activated in-cloud aerosol mass mixing ratio to domain top';                         
     case 'Total aerosol mass no act column integrated z1500m budget'
        ylab='Domain-mean (kg m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated aerosol mass without activated aerosol mixing ratio to z=1500m';                 
    case 'Coarse mode aerosol number column integrated z1500m budget'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated coarse mode aerosol number concentration to z=1500m';                   
    case 'Accum mode aerosol number column integrated z1500m budget'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated accumulation mode aerosol number concentration to z=1500m';                   
    case 'Total aerosol number column integrated budget z1500m to domain top'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated aerosol number concentration z=1500m to domain top';                   
    case 'droplet_number_column_integrated'
        ylab='Domain-mean (m^{-2})';
        titlenam='Timeseries of domain-mean of column integrated droplet number concentration';           
    otherwise
        error('Set titlenam, etc. here!');
end

if i_restricted_domain==1
    titlenam = [titlenam ' (restricted domain)'];
end

if iplot_flux_contribution==1
    ylab = [ylab ' FLUX contribution'];
end

% Shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 12th Nov this was at 16:48 for
% -20, -76 lat lon (centre of the domain). I.e. they are 4hrs 48 mins behind UTC
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)
xlims=1;
xlimits=[datenum('00:00 12-Nov-2008') datenum('00:00 16-Nov-2008')] + time_shift; %shift to LST since final plot will be in LST
xlimits=[datenum('00:00 12-Nov-2008') datenum('00:00 14-Nov-2008')] + time_shift; %shift to LST since final plot will be in LST

izlim=0;
zmin=1500;
zmax=3000;

lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=6; %=6; -1 means that all points have markers. Otherwise only plot the number specified.



idate_ticks_fix=1;
iaxis_square=0; %switch to make axis square

idat_micro=[];


%% -------- UM data -----------------------------------------------------

switch day_or_night
    case 'day'
        time = [ datenum('13-Nov-2008 17:00')  datenum('13-Nov-2008 23:00') ];
    case 'night'
        time = [ datenum('13-Nov-2008 02:00')  datenum('13-Nov-2008 08:00') ];
    case 'all'
        time = [ datenum('12-Nov-2008 06:00')  datenum('14-Nov-2008 00:00') ];
end

%% Script to get the UM run details by providing the run set name
%% Provide the case in UM_case_select_runs
UM_case_select_RUN  %runs UM_case_select_runs


clear LWP_overall_mean
        

%Bypass this for now
for idat_UM=1:length(fileUM)
%for idat_UM=1:length(fileUM)


    
    if iscell(dirUM)==1
        dirUM_i = dirUM{idat_UM};
    else
        dirUM_i = dirUM;
    end
    
    
    filename = [dirUM_i fileUM{idat_UM}];
    filename_rho = [dirUM_i fileUM_rho{idat_UM}];
    


%    lwp = get_LWP_RWP_UM('LWP',flag{idat_UM},filename,filename_rho);
%     clear vars_in
%     vars_in.var = 'LWP';
%     vars_in.flag = flag{idat_UM};
%     vars_in.file_lwp = filename;
%     vars_in.file_rho = filename_rho;
%     vars_in.pole_lat = pole_lat;
%     vars_in.pole_lon = pole_lon;
%     vars_in.time_in = [];
% 
%     [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
%     
    
%          vars_in.var = 'accum_number_total_column_to_z3000';
%          vars_in.flag = 'load_mat';
%          vars_in.file_lwp =  [dirUM_i remove_character(fileUM{idat_UM},'VAR_NAME','accum_number_total_column_to_z3000') '.mat'];   %Keep as file_lwp for loading .mat files
% %         vars_in.file_rho = [dirUM_i fileUM_rho{idat_UM}]; %filename_rho;
%          vars_in.pole_lat = pole_lat;
%          vars_in.pole_lon = pole_lon;
%          vars_in.time_in = []; %set to this for all times
         
%     [dat_UM2D_load,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);





switch var_tim
    case 'Total aerosol number column integrated budget'
        ydat_import(idat_driver).y = 0;
        VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        VAR_NAMES = {'accum_number','droplet_number'};       
        
        UM_add_timeseries
        
    case 'Total aerosol number column integrated z3000m budget'
        ydat_import(idat_driver).y = 0;
        %VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        %VAR_NAMES = {'accum_number','droplet_number'};       
        VAR_NAMES = {'accum_number_total_column_to_z3000','droplet_number_total_column_to_z3000','coarse_number_total_column_to_z3000'};
        
        UM_add_timeseries
        
    case 'Coarse mode aerosol number column integrated z1500m budget'
        ydat_import(idat_driver).y = 0;
        %VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        %VAR_NAMES = {'accum_number','droplet_number'};
        VAR_NAMES = {'coarse_number_total_column_to_z1500'};
        
        UM_add_timeseries
                
            
        
    case 'Total aerosol number column integrated z1500m budget'
        ydat_import(idat_driver).y = 0;
        %VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        %VAR_NAMES = {'accum_number','droplet_number'};
        
        if i_aerosol_processing_multi{idat_UM}==1        
            VAR_NAMES = {'accum_number_total_column_to_z1500','droplet_number_total_column_to_z1500','coarse_number_total_column_to_z1500'};
        else
            VAR_NAMES = {'accum_number_total_column_to_z1500','coarse_number_total_column_to_z1500'};
        end
        
        UM_add_timeseries
        
     case 'Accum mode aerosol number column integrated z1500m budget'
        ydat_import(idat_driver).y = 0;
        %VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        %VAR_NAMES = {'accum_number','droplet_number'};
        VAR_NAMES = {'accum_number_total_column_to_z1500'};
        
        UM_add_timeseries
        
    case 'Total aerosol mass column integrated z1500m budget'
        ydat_import(idat_driver).y = 0;
        %VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        %VAR_NAMES = {'accum_number','droplet_number'};
        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500'};        

        UM_add_timeseries
        
      case 'Total aerosol mass column integrated to specified z'
        ydat_import(idat_driver).y = 0;
%        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['accum_mass_total_column_to_z' hgt_str],['coarse_mass_total_column_to_z' hgt_str],['act_mass_liq_total_column_to_z' hgt_str]};        

        UM_add_timeseries  
        
    case 'Air mass column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        VAR_NAMES = {['air_mass_total_column_to_z' hgt_str]};

        UM_add_timeseries
        
      case 'Accumulation mode aerosol mass column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['accum_mass_total_column_to_z' hgt_str]};

        UM_add_timeseries  
        
      case 'Aitken mode aerosol mass column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['aitken_mass_total_column_to_z' hgt_str]};

        UM_add_timeseries       
        
       case 'Aitken mode aerosol number column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['aitken_number_total_column_to_z' hgt_str]};

        UM_add_timeseries        
        
        
      case 'Coarse mode aerosol mass column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['coarse_mass_total_column_to_z' hgt_str]};

        UM_add_timeseries          

      case 'Activated aerosol mass column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['act_mass_liq_total_column_to_z' hgt_str]};

        UM_add_timeseries      
              

    case 'Droplet number column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['droplet_number_total_column_to_z' hgt_str]};

        UM_add_timeseries
        
    case 'qL column integrated to specified z'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {['LWP_total_column_to_z' hgt_str]};

        UM_add_timeseries
        
    case 'Total aerosol mass column integrated to domain top budget'    
        ydat_import(idat_driver).y = 0;
%        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'accum_mass_total_column_to_zdomain_top','coarse_mass_total_column_to_zdomain_top','act_mass_liq_total_column_to_zdomain_top'};        
%        VAR_NAMES = {'accum_mass_total_column_to_zdomain_top','act_mass_liq_total_column_to_zdomain_top'};                

        UM_add_timeseries
        
     case 'Total aerosol number column integrated to domain top budget'    
        ydat_import(idat_driver).y = 0;
        if i_aerosol_processing_multi{idat_UM}==1 
            VAR_NAMES = {'accum_number_total_column_to_zdomain_top','coarse_number_total_column_to_zdomain_top','droplet_number_total_column_to_zdomain_top'};
        else %with processing off the aerosol number is not depleted, so don't need to take droplets into account.
            VAR_NAMES = {'accum_number_total_column_to_zdomain_top','coarse_number_total_column_to_zdomain_top'};                        
        end
        UM_add_timeseries       
        
     case 'Coarse mode aerosol mass column integrated to domain top budget'    
        ydat_import(idat_driver).y = 0;
%        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'coarse_mass_total_column_to_zdomain_top'};        

        UM_add_timeseries   
 
    case 'Coarse mode aerosol number column integrated to domain top budget'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'coarse_number_total_column_to_zdomain_top'};

        UM_add_timeseries


    case 'Accum mode aerosol number column integrated to domain top budget'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'accum_number_total_column_to_zdomain_top'};

        UM_add_timeseries
         
           
    case 'Droplet number column integrated to domain top budget'
        ydat_import(idat_driver).y = 0;
        %        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'droplet_number_total_column_to_zdomain_top'};

        UM_add_timeseries
        
        
         
        
     case 'Accum mode aerosol mass column integrated to domain top budget'    
        ydat_import(idat_driver).y = 0;
%        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'accum_mass_total_column_to_zdomain_top'};        

        UM_add_timeseries 
        
          
     case 'Activated in-cloud aerosol mass column integrated to domain top budget'    
        ydat_import(idat_driver).y = 0;
%        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'act_mass_liq_total_column_to_zdomain_top'};        

        UM_add_timeseries       
  
        
 case 'Total aerosol mass no act column integrated z1500m budget'
        ydat_import(idat_driver).y = 0;
        %VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        %VAR_NAMES = {'accum_number','droplet_number'};
        %VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500'};
        VAR_NAMES = {'accum_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500'};        
        
        UM_add_timeseries
        
                
    case 'Total aerosol number column integrated budget z1500m to domain top'
        ydat_import(idat_driver).y = 0;
        %VAR_NAMES = {'accum_number','coarse_number','droplet_number','aitken_number'}; %
        %VAR_NAMES = {'accum_number','droplet_number'};
        VAR_NAMES = {'accum_number_total_column_z1500_to_top','coarse_number_total_column_z1500_to_top','droplet_number_total_column_z1500_to_top'};
        
        UM_add_timeseries
        
    otherwise
         ydat_import(idat_driver).y = 0;
        VAR_NAMES = {var_tim};
        
        UM_add_timeseries
        
        
%          file_tim = [dirUM_i remove_character(fileUM{idat_UM},'VAR_NAME',var_tim) '_timeseries.mat'];
%          [dat_UM_timser,flux_UM_timser] = load_UM_timser(file_tim,append_str_timser{idat_UM});
%          time_matlab = dat_UM_timser.time_UM;
%          % 1D data straight from mat file
%          ydat_import(idat_driver).y = dat_UM_timser.timeseries_UM;
end


     
     
%          vars_in.var = 'LWP';
%          vars_in.flag = 'load_mat';
%          vars_in.file_lwp =  [dirUM_i remove_character(fileUM{idat_UM},'VAR_NAME','LWP') '.mat'];   %Keep as file_lwp for loading .mat files
% %         vars_in.file_rho = [dirUM_i fileUM_rho{idat_UM}]; %filename_rho;
%          vars_in.pole_lat = pole_lat;
%          vars_in.pole_lon = pole_lon;
%          vars_in.time_in = []; %set to this for all times
%     
%      [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
% 
%      ilwp = find(lwp<thresh_LWP_DRIVER);
%      Nd(ilwp)=NaN;

% Stuff for restricting the region when using 2D fields
%     stime=size(Nd,1);
%     [iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges] = restrict_to_region_2D_lats(LAT_val_DRIVER,LON_val_DRIVER,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,stime);
%     Nd=permute(Nd,[2 3 1]);
%     a=NaN*ones(size(Nd));
%     a(iregion_lin)=0;
%     Nd=Nd+a;
%     Nd=permute(Nd,[3 1 2]);
    
% vars_in.var = 'RWP';
%     vars_in.flag = flag{idat_UM};
%     vars_in.file_lwp = filename;
%     vars_in.file_rho = filename_rho;
%     vars_in.pole_lat = pole_lat;
%     vars_in.pole_lon = pole_lon;
%     vars_in.time_in = [];

%     [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
% 
%     stime=size(rwp,1);
%     [iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges] = restrict_to_region_2D_lats(LAT_val_DRIVER,LON_val_DRIVER,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,stime);
%     rwp=permute(rwp,[2 3 1]);
%     a=NaN*ones(size(rwp));
%     a(iregion_lin)=0;
%     rwp=rwp+a;
%     rwp=permute(rwp,[3 1 2]); 
    
%    rwp=0;
    
%    tlwp = lwp+rwp;
 
%  if using 2D data
%    ydat_import(idat_driver).y = meanNoNan(Nd(:,:),2);
    

    
    
    %Calculations for mean overall LWP and Nd
%     it_overall = find(time_matlab >= time(1) & time_matlab <= time(end) );
%     
%     lwp_temp = meanNoNan(lwp(it_overall,:,:),1);
%     LWP_overall_mean(idat_driver) = meanNoNan(lwp_temp(:),1);    
%     
%     Nd_temp = meanNoNan(Nd(it_overall,:,:),1);
%     Nd_overall_mean(idat_driver) = meanNoNan(Nd_temp(:),1);        

%    time=nc{'t'}(:);
%    t0_str=nc{'t'}.time_origin{1};
%    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%    xdat_import(idat_UM).x = datenum(t0_str2) + time;   


%    xdat_import(idat_driver).x = time_matlab; %
    xdat_import(idat_driver).x = time_matlab + time_shift; % 
    

    line_pattern(idat_driver)=line_patternUM(idat_UM);  line_colour(idat_driver)=line_colourUM(idat_UM); marker_style(idat_driver)=marker_styleUM(idat_UM);
    labs_import(idat_driver).l = labs_UM(idat_UM).l;
    

        
    if iplot_flux_increments==1
           line_pattern(idat_driver).p='-';
           
           line_pattern(idat_driver+1).p='--';
           labs_import(idat_driver+1).l = [labs_UM(idat_UM).l ' FLUX'];
           xdat_import(idat_driver+1).x = xdat_import(idat_driver).x(2:end); %use one less since the first difference shoudl be ignored since the data is zero at t=0
           ydat_import(idat_driver+1).y = ydat_import_incs(idat_driver).y(2:end); 
           
%Could change ylab here?
    end    
    
    
    idat_driver = idat_driver + idat_driver_inc;

end




%% ---  Main script to do plots and save
savedir = savedir_driver;
isave_plot=0; %DRIVER_lineplot_watervap can save the plot, but will save once have made adjustments after DRIVER_lineplot_watervap
%
ichoose_styles=1;
%  - Run watervap*
DRIVER_lineplot_watervap

for i=idat_micro   %length(h)-2
    set(h(i).h,'linestyle','none');
end

% set(h(1).h,'marker','*');
% set(h(1).h,'markersize',20);
% uistack(h(1).h,'top');

%Joint together the REMSS satellite values with a line.
% x_all=[]; y_all=[];
% for i=idat_micro  %1:length(xdat)-2
%    x_all = cat(2,x_all,xdat(i).x);
%    y_all = cat(2,y_all,ydat(i).y);   
% end
% [x_all,I]=sort(x_all);
% y_all=y_all(I);
% 
% plot(x_all,y_all,'b','linewidth',3);


switch var_tim 
    case 'accum_number_total_column_to_z3000';
%        set(gca,'ylim',[5e12 6e12]);
    case 'droplet_number_total_column_to_z3000';
%        set(gca,'ylim',[5e8 6e12]);

    case 'Total aerosol number column integrated budget'
       %set(gca,'ylim',[5e12 6e12]);
end

%Change the size of the window
pos=get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 1200 pos(4)]); 




% isave_plot=isave_plot_driver;
% if isave_plot==1
%     saveas_ps_fig_emf(gcf,[savename],'',0,1);
% end
