%% UKESM eval
% savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
% eval(['!mkdir ' savedir_date]); 

%Need to run 
%ACSIS_Robson_paper_load_data  - sets var_ukesm in there
%before this to load in the data. Also sets savedir_date. And sorts out the
%obs data (gcm_area_ etc.)

% var_ukesm = 'Nd_cf_weighted_UKESM';
% var_ukesm = 'calipso_low_cloud_amount';
% var_ukesm = 'SW_up_TOA';

regions={'1','3','8','4'};
regions={'4'};

iplot_inset=1;
iplot_individual_ens=0;


season='Annual';

%Can choose different periods for the trends
yr_start_trend_box = [1850 1970 2001];
yr_end_trend_box = [1970 2014 2014];
iplot_trend=[1 1 0]; %whether to plot the straight trend lines

%yr_start_trend_box = [1982]; %CCI CF
%yr_start_trend_box = [1984]; %Deep-C
%yr_start_trend_box = [2007];
%yr_end_trend_box = [2014];
%iplot_trend=[1]; %whether to plot the straight trend lines

ioverride_obs_trend_yr_start=0; %if =1 then uses yr_trend_obs_start otherwise uses start of 
%obs period
yr_trend_obs_start = yr_start_trend_box;

yr_start_plot = 1982; %where to start the plot from
yr_start_plot = 1844; %where to start the plot from
yr_title = 1965; %Year to plot the title at.
%yr_title = (yr_start_plot + 2014) / 2;

iplot_lin = 1; %whether to plot straight lines for the above linear trends

loc_DRIVER='SouthWest';
loc_DRIVER='NorthEast';

clear stats_regions trend_dat_box
for ibox=1:length(regions)

%% Choose region
% Box for regional averages
%box_region_DRIVER = '1'; %US small box
%box_region_DRIVER = '3'; %SW of Spain
%box_region_DRIVER = '4'; %All NA region
%box_region_DRIVER = '5'; %UK region
%box_region_DRIVER = '6'; %West of UK region
%box_region_DRIVER = '7'; %20 deg further west of there
%box_region_DRIVER = '8'; %middle of Atlantic at US box latitude
%box_region_DRIVER = '9'; %20 deg further west of there

box_region_DRIVER = regions{ibox};
%ACSIS_Robson_paper_plot_timeseries
ACSIS_Robson_paper_plot_timeseries3 %More streamlined modular version (Mar 2020).

end

%iscreen_land=0;

%%

