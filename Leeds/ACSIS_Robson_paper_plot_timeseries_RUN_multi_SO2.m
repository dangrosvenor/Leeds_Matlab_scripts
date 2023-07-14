%% UKESM eval
% savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
% eval(['!mkdir ' savedir_date]); 

%Need to run 
%ACSIS_Robson_paper_load_data  - sets var_ukesm in there
%before this to load in the data. Also sets savedir_date. And sorts out the
%obs data (gcm_area_ etc.)

season='Annual';
%season='DJF';
%season='MAM';
%season='JJA';
%season='SON';

%% Choose region
% Box for regional averages
%box_region_DRIVER = '1'; %'US outflow'  US small box
%box_region_DRIVER = '3'; %'Europe outflow' SW of Spain
%box_region_DRIVER = '4'; %'N. Atlantic basin' All NA region
%box_region_DRIVER = '5'; %UK region
%box_region_DRIVER = '6'; %West of UK region
%box_region_DRIVER = '7'; %20 deg further west of there
%box_region_DRIVER = '8'; %'Mid-Atlantic' middle of Atlantic at US box latitude
%box_region_DRIVER = '9'; %20 deg further west of there

regions={'1','3','8','4'};
regions={'4'};
%regions={'10'};
%regions={'11'};
regions={'12'};

iplot_individual_ens=0; %Set this if plotting the individual ensemble memebers
    % *** USE ACSIS_Robson_paper_plot_timeseries_RUN_multi_ENS_DeepC.m ***

%Can choose different periods for the trends
yr_start_trend_box = [1851 1970 2001];
yr_end_trend_box = [1970 2014 2014];
iplot_trend=[1 1 0]; %whether to plot the straight trend lines

yr_start_trend_box = [1985];
    switch season
        case 'DJF'
            yr_start_trend_box = yr_start_trend_box + 1; %start a year later for now - although could do DJF of Y2000 for CERES
    end

yr_end_trend_box = [2014];
iplot_trend=[1]; %whether to plot the straight trend lines

yr_start_plot = 1900; %where to start the plot from
yr_start_plot = 1850; %where to start the plot from
%yr_title = 1965; %Year to plot the title at.
yr_title = 1980; %Year to plot the title at.
%yr_title = (yr_start_plot + 2014) / 2;

yr_PI_bar = 1847.5;
yr_PI_bar = 2016;

iplot_lin = 1; %whether to plot straight lines for the above linear trends



clear stats_regions trend_dat_box
for ibox=1:length(regions)

box_region_DRIVER = regions{ibox};
%ACSIS_Robson_paper_plot_timeseries_deepc
ACSIS_Robson_paper_plot_timeseries_generic

end

%iscreen_land=0;

%%

