% -- Plot of SW vs Nd, along with the cloud properties (LWP, Nd and CF).
%  for the cloud properties will weight the averages by the mean SW_up, so
%  that the times of day that contribute most to the SW count more towards
%  the average.

%% mat file for saving the combined results from all the model runs
%% requested to (if isave_all_models==1)
isave_all_models = 0;
savefile_UM = ['/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/SW_values_for_ACI_' datestr(now,30) '.mat'];
savefile_UM = ['/home/disk/eos10/d.grosvenor/UM/Iceland/Nd_timeseries_' datestr(now,30) '.mat'];

%% Run for SW first - this will give us the weights for averaging

clear UM_time_in

%% Set the details on what data to get

%Choose the runs to plot/calculate for - see UM_case_select_runs
UM_time_in.UM_cases = '12th Nov case, as of May 2016';
%UM_time_in.UM_cases = '12th Nov case, as of May 2016 adhoc';
UM_time_in.UM_cases = 'Iceland_9day_runs_Nov2016';

%Choose the variable
%UM_time_in.varname = 'SW_up_TOA';
UM_time_in.varname = var_to_calc_multi{ivar_multi}; %from outside script
UM_time_in.var_units_str = ''; %may want to set this where set nc_varname, etc.

%Choose the lat lon required
UM_time_in.LAT_val_DRIVER = [-22.70 -17.28]; UM_time_in.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
UM_time_in.irestrict_domain_DRIVER=0; %set to zero to plot the full domain

%Choose times
  %For time specify either a range (will pick out all times within the range)
  %or a specific set of times (can be a single time). BUT DON'T SET NOT BOTH!!
  %UM_time_in.time_range =  [datenum('01-Sep-2014 12:00') datenum('01-Sep-2014 12:00')];
  %UM_time_in.time_specific(1) =  datenum('01-Sep-2014 12:00');
%UM_time_in.time_specific =  [datenum('12-Nov-2008 22:00') datenum('13-Nov-2008 12:00')];
UM_time_in.time_range =[ datenum('12-Nov-2008 06:00')  datenum('14-Nov-2008 00:00') ];
UM_time_in.time_range =[ datenum('12-Nov-2008 06:00')  datenum('14-Nov-2008 00:00') ];

%Iceland case
UM_time_in.time_range =  [datenum('31-Aug-2014 00:00') datenum('09-Sep-2014 21:02')];

%% Set the flags
UM_time_in.action='timeseries domain means';  %save_all
UM_time_in.iclose_figs=1; %whether to close the figures after plotting
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
[UM_out] = UM_maps_generic_time_loop_20141126T041216_v1(UM_time_in);



%% Save the results and the options used to get it
if isave_all_models == 1
    save(savefile_UM,'UM_out','UM_time_in','-V7.3');
end


