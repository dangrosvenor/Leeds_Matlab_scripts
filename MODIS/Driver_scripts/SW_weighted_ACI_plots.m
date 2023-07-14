% -- Plot of SW vs Nd, along with the cloud properties (LWP, Nd and CF).
%  for the cloud properties will weight the averages by the mean SW_up, so
%  that the times of day that contribute most to the SW count more towards
%  the average.

%% Run for SW first - this will give us the weights for averaging

clear UM_time_in

%% Set the details on what data to get

%Choose the runs to plot/calculate for 
UM_time_in.UM_cases = '12th Nov case, as of May 2016';

%Choose the variable
UM_time_in.varname = 'SW TOA out';
UM_time_in.var_units_str = 'W m^{-2}'; %may want to set this where set nc_varname, etc.

%Choose the lat lon required
UM_time_in.LAT_val_DRIVER = [-22.70 -17.28]; UM_time_in.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
UM_time_in.irestrict_domain_DRIVER=0; %set to zero to plot the full domain

%Choose times
  %For time specify either a range (will pick out all times within the range)
  %or a specific set of times (can be a single time). BUT DON'T SET NOT BOTH!!
  %UM_time_in.time_range =  [datenum('01-Sep-2014 12:00') datenum('01-Sep-2014 12:00')];
  %UM_time_in.time_specific(1) =  datenum('01-Sep-2014 12:00');
UM_time_in.time_specific =  [datenum('12-Nov-2008 22:00') datenum('13-Nov-2008 12:00')];


%% Set the flags
UM_time_in.isave_plot=0;
UM_time_in.iplot_maps=1;

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
[SW_UM] = UM_maps_generic_time_loop_20141126T041216_v1(UM_time_in);






