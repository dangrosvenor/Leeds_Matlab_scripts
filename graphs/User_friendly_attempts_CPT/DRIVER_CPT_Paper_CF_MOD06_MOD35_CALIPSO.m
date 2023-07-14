%creates the plots for the CF_MOD06_MOD35_CALIPSO comparison

%first need to load in AMSRE data using load_amsre_saved_data (select the
%2006-2010 Aqua only dataset).

%Then load in the MODIS data (Aqua only for the above years) using
%load_saved_modis_vars

%Then load the daytime CALIPSO CF data using
%read_calipso_monthly_daynight_IPSL

thresh_LAT = [-40 10];  thresh_LON = [-140 -50]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
thresh_CTH = 3.25; %km  - the approx height of 680 hPa (low cloud pressure threshold for ISCCP, CALIPSO, etc)
ifilter_clhcalipso=0; %don't do any cutting out of data based on the high CALIPSO CF
cal_CF_range='low'; %compare to low CALIPSO CF, or a combination?
cal_CF_range='low+mid+high';


%MOD35 CF without CTH screening
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot= 'Cloud_Fraction_Day_Mean'; mod_data_type='timeseries3 lambert'; %MOD35 daytime cloud fraction
%times_required = [12:15];  %Aqua daytime
screen_type ='none';  %screen_type ='CTH only';
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);

%MOD35 CF with CTH screening
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot= 'Cloud_Fraction_Day_Mean'; mod_data_type='timeseries3 lambert'; %MOD35 daytime cloud fraction
%times_required = [12:15];  %Aqua daytime
screen_type ='CTH only';
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);


%MOD06 without CTH screening
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot='Liquid Cloud Fraction time mean - specific days'; mod_data_type='timeseries3 lambert';
%times_required = [12:15];  %Aqua daytime
screen_type ='none';  %screen_type ='CTH only';
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);

%MOD06 for CTH < 3.25km only
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot='Liquid Cloud Fraction time mean - specific days'; mod_data_type='timeseries3 lambert';
%times_required = [12:15];  %Aqua daytime
screen_type ='CTH only';
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);


%Biases relative to CALISPO (CTH < 3.25km only).
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot='MODIS Cloud Fraction minus CALIPSO low CF'; mod_data_type='timeseries3 lambert';
%times_required = [12:15];  %Aqua daytime
screen_type ='CTH only';
MODCF_type = 'MOD35';
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);

ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot='MODIS Cloud Fraction minus CALIPSO low CF'; mod_data_type='timeseries3 lambert';
%times_required = [12:15];  %Aqua daytime
screen_type ='CTH only';
MODCF_type = 'MOD06 liquid';
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);


ifilter_clhcalipso=1;
thresh_clh = [0.3 1.01];
thresh_clh = [-0.01 0.3];

%Biases relative to CALISPO (CTH < 3.25km only).
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot='MODIS Cloud Fraction minus CALIPSO low CF'; mod_data_type='timeseries3 lambert';
%times_required = [12:15];  %Aqua daytime
screen_type ='CTH only';
MODCF_type = 'MOD35';
plot_global_maps
mean_bias_MOD35 = meanNoNan(P(:),1);
RMSE_MOD35 = sqrt(meanNoNan(P(:).^2,1));
saveas_ps_fig_emf(gcf,[savename],'',0,1);

ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
ioverride_plotglobal_loc=1;
time_mean_str ='ALL';
modis_data_plot='MODIS Cloud Fraction minus CALIPSO low CF'; mod_data_type='timeseries3 lambert';
%times_required = [12:15];  %Aqua daytime
screen_type ='CTH only';
MODCF_type = 'MOD06 liquid';
plot_global_maps
mean_bias_MOD06 = meanNoNan(P(:),1);
RMSE_MOD06 = sqrt(meanNoNan(P(:).^2,1));
saveas_ps_fig_emf(gcf,[savename],'',0,1);


