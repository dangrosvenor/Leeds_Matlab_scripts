% Gathers data on trends and for emission contributions.:
% Puts all the values into table_vals (loads from the .mat files)

% Should be able to run this script directly once all the trends have been
% calculated (i.e., without having to redo them)
% This script runs :-
% ACSIS_Robson_paper_TABLE_stats_noobs3_FUNC - may need to put entries in
%      here to read new variables or new runs.

%If do need to recalculate them :-
% run this script :-
% ACSIS_Robson_paper_load_data_plot_MULTI.m
% The above sets obs_str='none' before running the ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic.m script for each variable
% to save the data needed for this script.



%% Gather all the trend values and put into .tex file for table in the paper
clear table_vals table_vals_region0

idelta_vals=1; %whether to output the delta values (trends * delta time) or the trends


land_ocean_str = 'ocean_only';
land_ocean_str = 'ocean_only,_no_sea-ice';

ens_str=''; %previously had named them ensemble_mean, but seem to have stopped doing that.
iobs=0;

% -- Period 1
itr=1;
iobs_period=0;
period_str='PA';

%Run script for PA
ACSIS_Robson_paper_TABLE_stats_noobs3_FUNC

%% -- Period 2 --
itr=2; %specifies the second trend in the cell
iobs_period=0;
%i0=find(trends.years_ukesm_1d==trends.yr_start_trend_box2(itr));
period_str='PB';
%deltaT = trends.yr_end_trend_box2(itr) - trends.yr_start_trend_box2(itr);

%Run script for PB
ACSIS_Robson_paper_TABLE_stats_noobs3_FUNC

%% -- Period 3 --
% Some of the files are missing for the 3rd period (obs period - e.g.,
% HADGEM, DAMIP ones for Nd at least).
% So will just extract what is needed manually for now.

idelta_vals=0; %whether to output the delta values (trends * delta time) or the trends
    %Set to use the trends for the obs period since the delta T varies for
    %the different variables
itr=3; %specifies the second trend in the cell
iobs_period=1;
%i0=find(trends.years_ukesm_1d==trends.yr_start_trend_box2(itr));
period_str='PC';
%deltaT = trends.yr_end_trend_box2(itr) - trends.yr_start_trend_box2(itr);

%Run script for PC
ACSIS_Robson_paper_TABLE_stats_noobs3_FUNC

% %Want the UKESM and HADGEM model temperature deltas for the DEEP-C period
% %And the observed dT.
% %Then can calculate the difference in SW when use the observed dT compared
% %to the model ones (and later plot on the obs bar and whisker plots).
% %Will have the SW values there, so don't need to output those.
% %Will also need to calculate the dSW/dT as done for the bar plots.
% run_str='';
% var_str_mat='Surface_Temperature'; 
% var_str_tab='TS';
% deltaT=1; %Don't multiply by the time delta since the obs periods vary between variables.
% 
% hist_str='';
% tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
% trends = load(tr_file);
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
% 
% iobs=1;
% run_str='OBS';
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
% iobs=0;
% 
% %HADGEM
% run_str='HAD';
% hist_str='';
% tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
% trends = load(tr_file);
% ido_ens=1;
% 
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
% 
% tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
% %Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
% trends = load(tr_file);
% ido_ens=0;
% 
% %Hist-GHG
% run_str=['HistGhg'];
% hist_str='_hist_GHG';
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
% 
% 


% %Run for SW for DAMIP-hist-GHG since this calculates dSW/dT
% var_str_mat='SW_TOA_up'; 
% var_str_mat='F_SW_TOA_upwelling'; 
% var_str_tab='SW';
% 
% % -- DAMIP trends and uncertainties --
% %tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
% %tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
% tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
% %Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
% trends = load(tr_file);
% ido_ens=0;
% 
% %Hist-GHG
% run_str=['HistGhg'];
% hist_str='_hist_GHG';
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
% 
% 
% run_str='';
% hist_str='';
% tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
% trends = load(tr_file);
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
% 
% iobs=1;
% run_str='OBS';
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
% iobs=0;



%% save to .tex file

iappend=0;
save_sw_table_vals = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Trends_noobs_for_region_4_' land_ocean_str '_TABLE'];
var_str = '';
latex_newcommand_from_structure(table_vals,var_str,save_sw_table_vals,iappend);


