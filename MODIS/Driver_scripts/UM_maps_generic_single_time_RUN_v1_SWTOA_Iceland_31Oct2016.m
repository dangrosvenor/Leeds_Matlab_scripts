savefile = '/home/disk/eos10/d.grosvenor/UM/Iceland/saved_data_UM_maps_generic_time_loop_RUN_v1_SWTOA_Iceland_31Oct2016.mat';
isave_UM=0;

%options that will set for individual cases
%clear UM_time_in

%% Set the details on what data to get

% Set the action required (e.g. output time mean of all 2D fields).
% If not set then does 'output all', which just outputs the full 2D fields
% for each time requested.
%UM_time_in.UM_time_loop_action='time_mean_2D';
UM_time_in.UM_time_loop_action='output all';

%Choose the runs to plot/calculate for - see UM_case_select_runs
UM_time_in.UM_cases = '12th Nov case, as of May 2016';
UM_time_in.UM_cases = 'Iceland_9day_runs_Nov2016';
UM_time_in.UM_cases = 'Iceland_9day_runs_Nov2016_low_background_only';

%Choose the variable
%UM_time_in.varname = 'LWP';
UM_time_in.varname = 'SW_up_TOA';
UM_time_in.varname = 'SW_TOA_outgoing';
UM_time_in.var_units_str = 'W m^{-2}'; %may want to set this where set nc_varname, etc.

%Choose the lat lon required
UM_time_in.LAT_val_DRIVER = [-22.70 -17.28]; UM_time_in.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
UM_time_in.irestrict_domain_DRIVER=0; %set to zero to plot the full domain

% For shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 12th Nov this was at 16:48 for
% -20, -76 lat lon (centre of the domain). I.e. they are 4hrs 48 mins behind UTC
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)
  %N.B. not used below, just calculating here for potential future use

%% Choose times
% For time specify either a range (will pick out all times within the range)
% or a specific set of times (can be a single time). BUT DON'T SET NOT BOTH!!
% UM_time_in.time_range =  [datenum('01-Sep-2014 12:00') datenum('01-Sep-2014 12:00')];
% UM_time_in.time_specific(1) =  datenum('01-Sep-2014 12:00');
  
%UM_time_in.time_specific =  [datenum('13-Nov-2008 07:00')];  %For 02:12 LST snapshots for nightime figures
%UM_time_in.time_specific =  [datenum('12-Nov-2008 22:00') datenum('13-Nov-2008 12:00')];

%UM_time_in.time_specific =  [datenum('01-Sep-2014 12:01')];  %For 02:12 LST snapshots for nightime figures
%Times all seem to be one minute past hour...?
%UM_time_in.time_range =  [datenum('31-Aug-2014 00:00') datenum('09-Sep-2014 21:02')];

%UM_time_in.time_tol = 2/60/24; %set to 2 mins since times all seem to be 1 minute off


%% Set the flags
UM_time_in.isave_plot=0;
UM_time_in.iplot_maps=1;
UM_time_in.iclose_figs=1;

UM_time_in.i_mask_low_LWP=0; %Make any values below thresh_LWP equal to NaN
UM_time_in.thresh_LWP_mask = 20;

%Stuff for coarsening
UM_time_in.icoarsen = 0;
%UM_time_in.dlat_target=dlat_target; UM_time_in.dlon_target = dlon_target;

%Plotting options    
UM_time_in.iset_min_clim=1;
UM_time_in.clim_min=0;
UM_time_in.iset_max_clim=1;
UM_time_in.clim_max=650;

UM_time_in.iplot_markers=0; %marker points on the map

%% Run the time looping function
[UM_time_out] = UM_maps_generic_time_loop_20141126T041216_v1(UM_time_in);

if isave_UM==1
    save(savefile,'UM_time_out','UM_time_in');
end



