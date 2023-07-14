figlab='SW_trends_4panels';

%fig_dir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20200921T085951/';
fig_dir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';


clear fig_load pos_sub pos_cb

% --------------------------------------2D_PDF---
% Run this to set the defaults 
  defaults_multi_panel_plot
% -----------------------------------------

%N.B. Need to set color limits before this otherwise colobar will not
%change - specify in the script with the figure names in
Nsub = 2;
Msub = 2;  %No. of panels Nsub = no. vert, Msub = no. horiz

clims01 = [-25 25];

%stretch_cb=1; %whether to stretch the CB across the bottom of the plot
%fcb_stretch = 0.25; %1.7; %2.2

%Note for size_fac doesn't seem to do anything until get below a certain
%size...
size_fac = 0.93; %0.65; %0.85;
size_fac = 1;
%dY = [0.05 0.05 0.11 0.11]; %-0.02;  %0.05 %distance to move the subplot upwards relative to papersize
%dX = [0 -0.05 0 -0.05]; %0.4175;
dY = [0 0 0.05 0.05]; %-0.02;  %0.05 %distance to move the subplot upwards relative to papersize
dX = [0 -0.05 0 -0.05]; %0.4175;
% extra_for_even = 0; %flag to tell it to only move the odd numbered plot by dX
% dX_even = -0.1; %move the even plots by this much extra

%Not using these anymore
%dY_cb = 0.11; %0.01;  %0.11; %distance to move the colorbars upwards
%dX_cb = 0.0; %0.05

%dY_spacing = 0.17; %0.02; %vertical spacing to add between plots, negative values add space

%move the a,b,c labels in x and y by:-
dX_labs = -0.02;
dX_labs = 0.02;
dX_labs = -0.04; %-0.02; %-0.05;

dY_labs = -0.1;
dY_labs = -0.05;
dY_labs = 0.01; %0.07; %0.02; 


% ----  ------
ifig=1;
clims{ifig} = clims01;
pos_cb_set{ifig}=[0.114 0.1800 0.3500 0.0400];
fig_tit{ifig} = 'UKESM1 1850-1970';
fig_load{ifig} = 'global_UKESM1_rsut_trend_of_ensemble_mean_between_1850_and_1970;_W_m^-2.fig';
idelete_cb{ifig}=1;
idelete_xlab{ifig}=1;
idelete_xtick_labs{ifig}=1;
ioverride_xticks{ifig}=[1:8];
ioverride_yticks{ifig}=[1:8];
ifig=ifig+1;

clims{ifig} = clims01;
pos_cb_set{ifig}=[0.5520 0.1800 0.3500 0.0400];
fig_tit{ifig} = 'UKESM1 1971-2014';
fig_load{ifig} = 'global_UKESM1_rsut_trend_of_ensemble_mean_between_1971_and_2014;_W_m^-2.fig';
idelete_cb{ifig}=1;
idelete_xlab{ifig}=1;
idelete_ylab{ifig}=1;
idelete_xtick_labs{ifig}=1;
ioverride_xticks{ifig}=[1:8];
ioverride_yticks{ifig}=[1:8];
ifig=ifig+1;

clims{ifig} = clims01;
pos_cb_set{ifig}=[0.114 0.1800 0.3500 0.0400];
fig_tit{ifig} = 'HADGEM 1850-1970';
fig_load{ifig} = 'global_HADGEM-LL_rsut_trend_of_ensemble_mean_between_1850_and_1970;_W_m^-2.fig';
idelete_cb{ifig}=1;
ioverride_xticks{ifig}=[1:8];
ioverride_yticks{ifig}=[1:8];
ifig=ifig+1;

clims{ifig} = clims01;
pos_cb_set{ifig}=[0.150 0.09 0.690 0.0400];
fig_tit{ifig} = 'HADGEM 1971-2014';
fig_load{ifig} = 'global_HADGEM-LL_rsut_trend_of_ensemble_mean_between_1971_and_2014;_W_m^-2.fig';
%idelete_cb{ifig}=1;
idelete_ylab{ifig}=1;
ioverride_xticks{ifig}=[1:8];
ioverride_yticks{ifig}=[1:8];
icbar_horiz{ifig}=1; %Make the colorbar a horizontal one rather than vertical
cbar_label{ifig} = '\DeltaF_{\uparrowSW} (W m^{-2})';
ifig=ifig+1;





%Think can't make this name too long, or it will fail (Matlab limit)
run_multi_file ='run_multi_2D_PDF_vs_states';
eval(run_multi_file);
%finally, mnake a copy of the above file with the tag to match the plots
%(from saveas...). This will allow the figure to re-created in case
%run_multi... changes at any point. E.g. for different layouts. Will need
%to copy the logged file over the original run_multi* file (after backing it up) 
%eval_str = ['!C:\cygwin\bin\cp "C:\Users\Dan\Documents\MATLAB\work\graphs\' run_multi_file '.m" "C:\Users\Dan\Documents\MATLAB\work\graphs\run_multi_panel_copy_commands_log\' run_multi_file '_' datestr_now '.m"'];
%eval(eval_str);