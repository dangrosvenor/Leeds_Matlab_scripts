figlab='box_whisker_trends';

sim_period='Original'; %2009-2010 as used in ACPD verison of paper

fig_dir = '/home/disk/eos5/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';



clear fig_load pos_sub pos_cb

% --------------------------------------2D_PDF---
% Run this to set the defaults 
  defaults_multi_panel_plot
% -----------------------------------------

%N.B. Need to set color limits before this otherwise colobar will not
%change - specify in the script with the figure names in
Nsub = 3;
Msub = 2;  %No. of panels Nsub = no. vert, Msub = no. horiz

%stretch_cb=1; %whether to stretch the CB across the bottom of the plot
%fcb_stretch = 0.25; %1.7; %2.2

%Note for size_fac doesn't seem to do anything until get below a certain
%size...
size_fac_X = 0.93; %0.65; %0.85;
size_fac_Y = 0.83; %0.65; %0.85;

%dX = [0 0 0 0 0 0 0 0]; 
%dY = [0 0 0 0 0 0 0 0];

%dX = [0.03 0.03 0.03 0.03 0.03 0 0 0]+0.03; %-0.02;  %0.05 %distance to move the subplot upwards relative to papersize
dX = [0.035 0.08 0.035 0.075 0.04 0.075 0 0]+0.03; %-0.02;  %0.05 %distance to move the subplot upwards relative to papersize
dY = [0.03 0.03 0.005 0.005 -0.02 -0.02 0 0]; %-0.02;  %0.05 %distance to move the subplot upwards relative to papersize

%dY = [0.05 0.05 0.11 0.11]; %-0.02;  %0.05 %distance to move the subplot upwards relative to papersize
%dX = [0 -0.05 0 -0.05]; %0.4175;


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
dY_labs = 0.0; %0.07; %0.02; 

%Adjust position of xlabel
dX_xlab=[0 0 0 0 0 0];
dY_xlab=[0.15 0.18 0.10 0.15 0.15 0];
dY_xlab=[0 0 0 0 0 0];
%dY_xlab=[0.12 0.12 0.12 0.12 0.12 0];
%dY_xlab=0.0; %0.4
%fsize_xlab = 16;



ifig=1;
clims{ifig} = [-2 2];
pos_cb_set{ifig}=[0.114 0.1800 0.3500 0.0400];
fig_tit{ifig} = '';
%fig_load{ifig} = 'Annual_mean_SW_TOA_up_for_region_4_ensemble_mean,__ocean_only_box_whisk.fig'; 
%fig_load{ifig} = 'Annual_mean_F_SW_TOA_upwelling_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig'; 
fig_load{ifig} = 'Annual_mean_F_SW_uparrow_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
idelete_cb{ifig}=1;
idelete_xlab{ifig}=0;
idelete_xtick_labs{ifig}=0;
%ioverride_xticks{ifig}=[1:8];
%ioverride_yticks{ifig}=[1:8];
ifig=ifig+1;

clims{ifig} = [-2 2];
pos_cb_set{ifig}=[0.114 0.1800 0.3500 0.0400];
fig_tit{ifig} = '';
%fig_load{ifig} = 'Annual_mean_N_d_for_region_4_ensemble_mean_ocean_only_box_whisk.fig'; 
fig_load{ifig} = 'Annual_mean_N_d_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
idelete_cb{ifig}=1;
%ioverride_xticks{ifig}=[1:8];
%ioverride_yticks{ifig}=[1:8];
ifig=ifig+1;

clims{ifig} = [-2 2];
pos_cb_set{ifig}=[0.114 0.1800 0.3500 0.0400];
fig_tit{ifig} = '';
%fig_load{ifig} = 'Annual_mean_N_d_for_region_4_ensemble_mean_ocean_only_box_whisk.fig'; 
%fig_load{ifig} = 'Annual_mean_Aerosol+Dust_Optical_Depth,550nm_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
fig_load{ifig} = 'Annual_mean_tau_a_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
idelete_cb{ifig}=1;
%ioverride_xticks{ifig}=[1:8];
%ioverride_yticks{ifig}=[1:8];
ifig=ifig+1;

clims{ifig} = [-2 2];
pos_cb_set{ifig}=[0.5520 0.1800 0.3500 0.0400];
fig_tit{ifig} = '';
%fig_load{ifig} = 'Annual_mean_totCF_for_region_4_ensemble_mean,__ocean_only_box_whisk.fig'; 
%fig_load{ifig} = 'Annual_mean_totCF_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
%fig_load{ifig} = 'Annual_mean_totCF_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
fig_load{ifig} = 'Annual_mean_f_c_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
idelete_cb{ifig}=1;
idelete_xlab{ifig}=0;
idelete_ylab{ifig}=0;
idelete_xtick_labs{ifig}=0;
%ioverride_xticks{ifig}=[1:8];
%ioverride_yticks{ifig}=[1:8];
ifig=ifig+1;



clims{ifig} = [-2 2];
pos_cb_set{ifig}=[0.150 0.09 0.690 0.0400];
fig_tit{ifig} = '';
%fig_load{ifig} = 'Annual_mean_LWP_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig'; 
fig_load{ifig} = 'Annual_mean_L_all-sky_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
%idelete_cb{ifig}=1;
%idelete_ylab{ifig}=1;
%ioverride_xticks{ifig}=[1:8];
%ioverride_yticks{ifig}=[1:8];
%icbar_horiz{ifig}=0; %Make the colorbar a horizontal one rather than vertical
cbar_label{ifig} = '';
ifig=ifig+1;

% clims{ifig} = [-2 2];
% pos_cb_set{ifig}=[0.150 0.09 0.690 0.0400];
% fig_tit{ifig} = '';
% fig_load{ifig} = 'Annual_mean_LWPic_for_region_4_ensemble_mean,__ocean_only_box_whisk.fig'; 
% %idelete_cb{ifig}=1;
% %idelete_ylab{ifig}=1;
% %ioverride_xticks{ifig}=[1:8];
% %ioverride_yticks{ifig}=[1:8];
% %icbar_horiz{ifig}=1; %Make the colorbar a horizontal one rather than vertical
% %cbar_label{ifig} = '';
% ifig=ifig+1;

clims{ifig} = [-2 2];
pos_cb_set{ifig}=[0.150 0.09 0.690 0.0400];
fig_tit{ifig} = '';
%fig_load{ifig} = 'Annual_mean_Surface_Temperature_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig'; 
fig_load{ifig} = 'Annual_mean_T_for_region_4__ocean_only,_no_sea-ice_UKESM1__box_whisk.fig';
%idelete_cb{ifig}=1;
%idelete_ylab{ifig}=1;
%ioverride_xticks{ifig}=[1:8];
%ioverride_yticks{ifig}=[1:8];
%icbar_horiz{ifig}=1; %Make the colorbar a horizontal one rather than vertical
%cbar_label{ifig} = '';
ifig=ifig+1;



%Think can't make this name too long, or it will fail (Matlab limit)
run_multi_file ='run_SW_paper_trend_bars';
eval(run_multi_file);
%finally, mnake a copy of the above file with the tag to match the plots
%(from saveas...). This will allow the figure to re-created in case
%run_multi... changes at any point. E.g. for different layouts. Will need
%to copy the logged file over the original run_multi* file (after backing it up) 
%eval_str = ['!C:\cygwin\bin\cp "C:\Users\Dan\Documents\MATLAB\work\graphs\' run_multi_file '.m" "C:\Users\Dan\Documents\MATLAB\work\graphs\run_multi_panel_copy_commands_log\' run_multi_file '_' datestr_now '.m"'];
%eval(eval_str);