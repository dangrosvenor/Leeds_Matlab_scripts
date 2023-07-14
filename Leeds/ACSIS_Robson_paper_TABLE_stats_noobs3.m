% Should be able to run this script directly once all the trends have been
% calculated (i.e., without having to redo them)
% This script produces bars plots for multiple variables for the emission
% contributions (aeroso, GHG, etc.).
% It also does the feedback bar plots.

% This script runs :-
% ACSIS_Robson_paper_make_TABLE_emisison_contributions
%    which in turn runs
% ACSIS_Robson_paper_TABLE_stats_noobs3_FUNC - may need to put entries in
%    there to read new variables or new runs.
% ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot - plots the bar/errorbar
% plots - SET the variables to plot for here.

%If do need to recalculate them :-
% run this script :-
% ACSIS_Robson_paper_load_data_plot_MULTI.m
% The above sets obs_str='none' before running the ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic.m script for each variable
% to save the data needed for this script.

inew_folder=0;
if inew_folder==1
    savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
    eval(['!mkdir ' savedir_date]);
else
    savedir_date='/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';
end

%% Gather all the trend values and put into .tex file for table in the paper
ACSIS_Robson_paper_make_TABLE_emisison_contributions


%% Plot bar plots
ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot

