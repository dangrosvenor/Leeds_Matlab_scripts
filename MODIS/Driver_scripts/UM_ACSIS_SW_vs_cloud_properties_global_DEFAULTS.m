% N.B. - to get a global plot can just set irestrict_domain_DRIVER=0

%Pre-process using
%   UM_quick_plot_global.m

%Add variables to climits are set here too :-
%   UM_var_defs.m

%runs this plotting script - lat/lon for map also set here. FOR PDFs, etc. it is set below :-
  %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global  
  


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

%[out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,dat_global.time_ALL,time_choice,dim);

iload_rwp=1;
iload_SW_down_TOA=1;
amsre_data = 1;

ical_data=1; %Whether are using CALIPSO data - load using script :- read_calipso_monthly_night_IPSL.m for average values
    %first. This puts the data into these fields :- 
    % cllcalipso_monthly_AVERAGE, clmcalipso_monthly_AVERAGE, clhcalipso_monthly_AVERAGE
    % Is average of day and night values for now (have separate ones too).

iceres_data=1; %whether to load CERES data

cont_col_str_DRIVER='k';
ibias_contour=0;
icontour_DRIVER=0;
