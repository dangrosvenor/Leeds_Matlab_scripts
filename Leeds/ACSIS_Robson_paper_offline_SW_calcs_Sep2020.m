%This one plots the estimated trend line for the full calculated trend and
%also those with the different vars held constant. Also saves the trend
%data from those in a .mat file that is used to make the table and bar charts (using
% ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE.m

% Will attempt to load a .mat file with the trends from the SW calc based on
% the years chosen and iconstant_trend - if it doesn't exist then need to run 
%   ACSIS_Robson_paper_offline_SW_calcs_Sep2020_CALCs.m
% to calculate and save the values.


inew_folder=0;
if inew_folder==1
    savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
    eval(['!mkdir ' savedir_date]);
else
    savedir_date='/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';
end


%% set models 
%MIP = 'CMIP'; DAMIP_runs = {''}; DAMIP_runs2 = {'ukesm'};
MIP = 'DAMIP'; DAMIP_runs = {'DAMIP_hist-aer','DAMIP_hist-GHG','DAMIP_hist-nat'}; DAMIP_runs2 = DAMIP_runs;
MIP = 'DAMIP'; DAMIP_runs = {'DAMIP_hist-GHG'}; DAMIP_runs2 = DAMIP_runs;
%MIP = 'DAMIP'; DAMIP_runs = {'DAMIP_hist-GHG'}; DAMIP_runs2 = DAMIP_runs;
MIP = 'AerChemMIP'; DAMIP_runs = {'AerChemMIP_hist-piAer'}; DAMIP_runs2 = DAMIP_runs;
%MIP = 'AerChemMIP'; DAMIP_runs = {'AerChemMIP_hist-AerProxy'}; DAMIP_runs2 = DAMIP_runs;
%MIP = 'AerChemMIP'; DAMIP_runs = {'UKESM1-AerChemMIP_control'}; DAMIP_runs2 = DAMIP_runs;
%MIP = 'HADGEM3_GC31_LL'; DAMIP_runs = {'HADGEM3_GC31_LL'}; DAMIP_runs2 = DAMIP_runs;

%% set period, etc.
%iconstant_trend=1; iplot_trend=[1 0]; yr_start_plot = 1850; yr_end_plot = 1980; period_str='_P1';
%iconstant_trend=2; iplot_trend=[0 1]; yr_start_plot = 1970;  yr_end_plot = 2018; period_str='_P2';

%iconstant_trend=1; iplot_trend=[1 1]; yr_start_plot = 1850;

yr_start_trend_box = [1870 1985]; yr_end_trend_box = [1970 2014]; %
yr_start_trend_box = [1850 1971]; yr_end_trend_box = [1970 2014]; %

%% Load the mask for seaice regions based on max annual surface albedo for DAMIP runs over all ens members
% Calculated using ACSIS_Robson_paper_Max_surf_albedo.m

max_surf_dir = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/'];
max_surf_name = [max_surf_dir 'max_surf_albedo_DAMIP_annual'];
load(max_surf_name);

%Next need to figure out how screen based on this data... prob do at the
%same time as landmask - could combine this mask to that one I guess

for iperiod=1:2
    
    if iperiod==1        
        iconstant_trend=1; iplot_trend=[1 0]; yr_start_plot = 1850; yr_end_plot = 1980; period_str='_P1';
    else        
        iconstant_trend=2; iplot_trend=[0 1]; yr_start_plot = 1979; yr_end_plot = 2018; period_str='_P2';               
    end
    
    
% Loop over various models


%Loop over all models

for idamip_run=1:length(DAMIP_runs)
        
    expt_str = DAMIP_runs{idamip_run};    
    expt_str2 = DAMIP_runs2{idamip_run};   

%% Load the calculated data from the .mat file
   
% savedir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';
% savename=[savedir 'offline_SW_calcs_start_' num2str(yr_start_trend_box) '_end_' num2str(yr_end_trend_box) '_trend_period_' ...
%     num2str(iconstant_trend) expt_str '.mat'];
% savename = remove_character(savename,' ','_');


%savedir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';

savedir = '/home/disk/eos15/d.grosvenor/UM/ACSIS_Nd_trends/';
load_mat_file=[savedir 'offline_SW_calcs_time_varying_start_' num2str(yr_start_trend_box) '_end_' num2str(yr_end_trend_box) '_trend_period_' ...
    num2str(iconstant_trend) expt_str '_post_calc_ens_mean.mat'];
load_mat_file = remove_character(load_mat_file,' ','_');

load(load_mat_file);




%/home/disk/eos15/d.grosvenor/UM/ACSIS_Nd_trends/offline_SW_calcs_time_varying_start_1850__1971_end_1970__2014_trend_period_1DAMIP_hist-aer_post_calc_ens_mean.mat

%% Create dat_ukesm for use in the timeseries plotting script - plot timeseries

%% Plot overall clear-sky SW - N.B. - estimated and observed will be the same if using the clear-sky fluxes directly
% in the calculations.

iplot_box_whisker=0;
f_sw_calc = 1.0; %default
switch MIP
    case 'CMIP'
        f_sw_calc=1.66;
        %f_sw_calc=1.73;
        f_sw_calc=1.0;
    case {'DAMIP','AerChemMIP'}
        f_sw_calc = 1.6956;
        f_sw_calc = 1.7;        
        %f_sw_calc = 0.8857;   
        %f_sw_calc = 0.89;  
        %f_sw_calc = 0.8879; %best match after changing trans to 0.83 - actually changing trans and this number are equivalent.
            %I think. Except that is propto trans.^2. So if wanted to
            %adjust trans would need to make it. sw1=f*sw0 =
            %(trans1/trans0).^2 *sw0  = (trans1/0.83).^2 *sw0. trans1 =
            %sqrt(f)*trans0 = sqrt(0.8879)*0.83 = 0.7821. Previously was
            %0.6
        %f_sw_calc = 0.8879/0.9834; %test to match the 2nd period better for Hist-Aer
        f_sw_calc = 1.0; 
        %f_sw_calc = 1.011; 
        %f_sw_calc = 1.011*0.5784*0.9812*1.005;
        %f_sw_calc = 0.8719;
        %f_sw_calc = 0.9267;
        %f_sw_calc = 0.935;
        %f_sw_calc = 1.21; 
        f_sw_calc = 92.65/95.13 * 91.85/98.8 * 94.57/92.15 * 96.23/89.46;
        
        f_sw_calc = 0.9093;
        f_sw_calc = 1.0;        
end


switch expt_str
    case '' %UKESM        
        model_str = 'UKESM1';
    otherwise
        model_str = expt_str; %expt_str='DAMIP_hist-aer', etc.
end

    

%Clear-sky TOA flux from model
var_DAMIP = 'rsutcs';
load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str2 '_all_' var_DAMIP '.mat'];
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm
mat_obj_cs = matfile(load_file);

% Now put in the data that loaded earlier
dat_ukesm.dat_annual = mat_obj_cs.dat_annual; %SW_up_TOA_dat_ens_mean;
dat_ukesm.dat_annual_ens = mat_obj_cs.dat_annual_ens; %Just used for the inter-ensemble spread

%Treat the calculated (offline) ones like the obs
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
dat_ukesm_swcalc = load(load_file,'gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_ukesm_swcalc.years_ukesm_1d = years_sw_calc';

dat_ukesm_swcalc.dat_annual = f_sw_calc * SWTOA_clear_sky_calc_model_annual;
dat_ukesm_swcalc.dat_monthly = f_sw_calc * SWTOA_clear_sky_calc_model; %are the monthly ones used?
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
% a = repmat(f_sw_calc * SWTOA_calc_model_annual,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm_swcalc.dat_annual_ens=a;


%Run the ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic to make the
%timeseries for the  SW calculation with all variables chagning for
%comparison to the actual UKESM SW.
%This will also save the trend data to the .mat file. period_str below gets
%added to the filename of the .mat file, so that we have different
%filenames for this, the P1 plot (using constant values from P1) and the P2
%plot.
iadd_trend_str=0;
iplot_box_whisker=0;
iremove_xlabs = 0;
ipad_legend = 0; %Pad the first entry to avoid box cutting into the text when have
                    %superscripts
iplot_trend=[1 1]; yr_start_plot = 1850; yr_end_plot = 2018; period_str=''; plot_year_lims = [1848 2018];
var_ukesm='SWTOA Calc';
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic


%% Plot overall SW

iplot_box_whisker=0;




switch MIP
    case 'CMIP'
        f_sw_calc=1.66;
        %f_sw_calc=1.73;
        f_sw_calc=1.0;
        f_sw_calc=0.977;
        
    case {'DAMIP','AerChemMIP'}
        f_sw_calc = 1.6956;
        f_sw_calc = 1.7;        
        %f_sw_calc = 0.8857;   
        %f_sw_calc = 0.89;  
        %f_sw_calc = 0.8879; %best match after changing trans to 0.83 - actually changing trans and this number are equivalent.
            %I think. Except that is propto trans.^2. So if wanted to
            %adjust trans would need to make it. sw1=f*sw0 =
            %(trans1/trans0).^2 *sw0  = (trans1/0.83).^2 *sw0. trans1 =
            %sqrt(f)*trans0 = sqrt(0.8879)*0.83 = 0.7821. Previously was
            %0.6
        %f_sw_calc = 0.8879/0.9834; %test to match the 2nd period better for Hist-Aer
        f_sw_calc = 1.0; 
        %f_sw_calc = 1.011; 
        %f_sw_calc = 1.011*0.5784*0.9812*1.005;
        %f_sw_calc = 0.8719;
        %f_sw_calc = 0.9267;
        %f_sw_calc = 0.935;
        %f_sw_calc = 1.21; 
        f_sw_calc = 92.65/95.13 * 91.85/98.8;
        f_sw_calc = 0.8638* 96.23/89.46;
        f_sw_calc = 0.9093*0.9577; %0.9139*0.995;
        f_sw_calc = 0.8612;
        f_sw_calc = 1;
        %f_sw_calc = 88.2/89.5;
    case 'HADGEM3_GC31_LL'
        f_sw_calc = 0.989;
        f_sw_calc = 0.98;
        f_sw_calc = 1;
end


switch expt_str
    case '' %UKESM        
        model_str = 'UKESM1';
    otherwise
        model_str = expt_str; %expt_str='DAMIP_hist-aer', etc.
end




    

%Put the SWTOA data from the UKESM (i.e., full online SW) in the dat_ukesm file
%First load in the lat/lon etc.
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

% Now put in the data that loaded earlier
dat_ukesm.dat_annual = SW_up_TOA_dat_annual_ens_mean; %SW_up_TOA_dat_ens_mean;
dat_ukesm.dat_annual_ens = SW_up_TOA_dat_annual_ens;

%test since actual SWTOA fluxes look weird
%mat_obj = matfile('/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_DAMIP_hist-aer_all_rsut.mat');
%dat_ukesm.dat_annual = mat_obj.dat_annual;
%dat_ukesm.dat_annual = meanNoNan(mat_obj.dat_annual_ens,1); %these two work out the same

%mat_obj2 = matfile(load_mat_file);
%dat_ukesm.dat_annual = mat_obj2.SW_up_TOA_dat_annual_ens_mean;

%could test teh individual ense members

%Treat the calculated (offline) ones like the obs
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
dat_ukesm_swcalc = load(load_file,'gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
dat_ukesm_swcalc.years_ukesm_1d = years_sw_calc';

dat_ukesm_swcalc.dat_annual = f_sw_calc * SWTOA_calc_model_annual;
dat_ukesm_swcalc.dat_monthly = f_sw_calc * SWTOA_calc_model;
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
% a = repmat(f_sw_calc * SWTOA_calc_model_annual,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm_swcalc.dat_annual_ens=a;


%Run the ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic to make the
%timeseries for the  SW calculation with all variables chagning for
%comparison to the actual UKESM SW.
%This will also save the trend data to the .mat file. period_str below gets
%added to the filename of the .mat file, so that we have different
%filenames for this, the P1 plot (using constant values from P1) and the P2
%plot.
iadd_trend_str=0;
iplot_box_whisker=0;
iremove_xlabs = 1;
ipad_legend = 0; %Pad the first entry to avoid box cutting into the text when have
                    %superscripts
iplot_trend=[1 1]; yr_start_plot = 1850; yr_end_plot = 2018; period_str=''; plot_year_lims = [1848 2018];
var_ukesm='SWTOA Calc';
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic



%% Choose the period to plot for the following (Nd, CF, LWPic)
itrend_lab = iconstant_trend;
switch iconstant_trend
    case 1
        iplot_trend=[1 0]; yr_start_plot = 1850; yr_end_plot = 1980; period_str='_P1';
        plot_year_lims = [1848 1972];
        %
    case 2
        iplot_trend=[0 1]; yr_start_plot = 1970;  yr_end_plot = 2018; period_str='_P2';
        plot_year_lims = [1970 2014];
end

%iadd_trend_str=1; %whether to add the trend value to the legend
iadd_trend_str=1;
iplot_box_whisker=0;

%% Nd constant (not used any more)

iremove_xlabs = 1;
ipad_legend = 1; %Pad the first entry to avoid box cutting into the text when have
                    %superscripts

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_Nd;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_Nd,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str = 'N_d constant';
    otherwise
        model_str = [expt_str '_N_d constant']; %expt_str='DAMIP_hist-aer', etc.
end
%model_str='N_d constant';
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic



%% CF constant (not used any more)
iremove_xlabs = 1;
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_cf;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_cf,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='CF constant';
    otherwise
        model_str = [expt_str '_CF constant']; %expt_str='DAMIP_hist-aer', etc.
end

obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic


%% LWP constant (not used any more)
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_lwp;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_lwp,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='LWP constant';
    otherwise
        model_str = [expt_str '_LWP constant']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic

%% Nd varying (all else constant)
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_Nd_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_Nd_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='Nd vary';
    otherwise
        model_str = [expt_str '_Nd vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic


%% LWPic varying (all else constant)
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_lwp_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_lwp_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='lwp vary';
    otherwise
        model_str = [expt_str '_lwp vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic


%% CF varying (all else constant)
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_cf_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_cf_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='cf vary';
    otherwise
        model_str = [expt_str '_cf vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic

%% Surface albedo varying (all else constant)
isurf_vary=1;
if isurf_vary==1
    
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_albedo_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_albedo_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='albedo vary';
    otherwise
        model_str = [expt_str '_albedo vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic

end


%% Transmissivity varying (all else constant)
itrans = 0;
if itrans==1

iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_trans_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_trans_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='trans vary';
    otherwise
        model_str = [expt_str '_trans vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic

end

%% AOD varying (all else constant)
iaod_vary = 0;
if iaod_vary==1
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_aod_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_aod_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='aod vary';
    otherwise
        model_str = [expt_str '_aod vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic

end

%% Absorbing AOD varying (all else constant)
iaaod_vary=0;
if iaaod_vary==1
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_aaod_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_aaod_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='aaod vary';
    otherwise
        model_str = [expt_str '_aaod vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic

end

%% Clear-sky flux varying (all else constant)
iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
ipad_legend = 0;

%Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%script to set up some things needexd in dat_ukesm, etc.
ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_clear_sky_vary;
%Make dummy data for the ensemble
a = repmat(f_sw_calc * SWTOA_calc_model_annual_clear_sky_vary,[1 1 1 2]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

var_ukesm='SWTOA Calc';
switch expt_str
    case '' %UKESM        
        model_str='clear sky vary';
    otherwise
        model_str = [expt_str '_clear sky vary']; %expt_str='DAMIP_hist-aer', etc.
end
obs_str='SW calc';
ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.

%ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic





%%

% %% cf and LWPic constant and vary Nd - using period mean CF and LWPic (Nd contribution).
% iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
% ipad_legend = 0;
% 
% %Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
% load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
%%script to set up some things needexd in dat_ukesm, etc.
%ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm
% 
% dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_cf_lwp;
% %Make dummy data for the ensemble
% a = repmat(f_sw_calc * SWTOA_calc_model_annual_cf_lwp,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm.dat_annual_ens=a;
% 
% var_ukesm='SWTOA Calc';
% 
% var_ukesm='SWTOA Calc';
% switch expt_str
%     case '' %UKESM        
%         model_str='cf_LWP constant';
%     otherwise
%         model_str = [expt_str '_cf_LWP constant']; %expt_str='DAMIP_hist-aer', etc.
% end
% obs_str='SW calc';
% ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.
% 
% %ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic
% 
% 
% %% cf and Nd constant and vary LWPic - using period mean CF and Nd (LWPic contribution).
% iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
% ipad_legend = 0;
% 
% %Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
% load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
% %script to set up some things needexd in dat_ukesm, etc.
% ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm
% 
% dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_Nd_cf;
% %Make dummy data for the ensemble
% a = repmat(f_sw_calc * SWTOA_calc_model_annual_Nd_cf,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm.dat_annual_ens=a;
% 
% var_ukesm='SWTOA Calc';
% 
% var_ukesm='SWTOA Calc';
% switch expt_str
%     case '' %UKESM        
%         model_str='Nd_cf constant';
%     otherwise
%         model_str = [expt_str '_Nd_cf constant']; %expt_str='DAMIP_hist-aer', etc.
% end
% obs_str='SW calc';
% ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.
% 
% %ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic
% 
% 
% %% Nd and LWPic constant and vary Nd - using period mean Nd and LWPic (cf contribution).
% iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
% ipad_legend = 0;
% 
% %Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
% load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
% %script to set up some things needexd in dat_ukesm, etc.
% ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm
% 
% dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_lwp_Nd;
% %Make dummy data for the ensemble
% a = repmat(f_sw_calc * SWTOA_calc_model_annual_lwp_Nd,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm.dat_annual_ens=a;
% 
% var_ukesm='SWTOA Calc';
% 
% var_ukesm='SWTOA Calc';
% switch expt_str
%     case '' %UKESM        
%         model_str='LWP_Nd constant';
%     otherwise
%         model_str = [expt_str '_LWP_Nd constant']; %expt_str='DAMIP_hist-aer', etc.
% end
% obs_str='SW calc';
% ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.
% 
% %ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic
% 
% 
% 
% 
% 
% %% cf and LWPic constant and vary Nd - CF and LWPic from the start of the trend periods
% iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
% ipad_legend = 0;
% 
% %Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
% load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
% %script to set up some things needexd in dat_ukesm, etc.
% ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm
% 
% dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_cf_lwp_Pstart;
% %Make dummy data for the ensemble
% a = repmat(f_sw_calc * SWTOA_calc_model_annual_cf_lwp_Pstart,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm.dat_annual_ens=a;
% 
% var_ukesm='SWTOA Calc';
% 
% var_ukesm='SWTOA Calc';
% switch expt_str
%     case '' %UKESM        
%         model_str='cf_LWP constant';
%     otherwise
%         model_str = [expt_str '_cf_LWP constant Pstart']; %expt_str='DAMIP_hist-aer', etc.
% end
% obs_str='SW calc';
% ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.
% 
% %ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic
% 
% 
% %% cf and LWPic constant and vary Nd - using CF and LWPic from end of periods
% iremove_xlabs = 0; %keep xlabs since is the last plot of the subfig
% ipad_legend = 0;
% 
% %Put the test (Nd constant, etc.) calculated data in the dat_ukesm file
% load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat'; load(load_file,'years_ukesm_1d'); 
% %script to set up some things needexd in dat_ukesm, etc.
% ACSIS_Robson_paper_offline_SW_calcs_Sep2020_setup_datukesm
% 
% dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual_cf_lwp_Pend;
% %Make dummy data for the ensemble
% a = repmat(f_sw_calc * SWTOA_calc_model_annual_cf_lwp_Pend,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm.dat_annual_ens=a;
% 
% var_ukesm='SWTOA Calc';
% 
% var_ukesm='SWTOA Calc';
% switch expt_str
%     case '' %UKESM        
%         model_str='cf_LWP constant';
%     otherwise
%         model_str = [expt_str '_cf_LWP constant Pend']; %expt_str='DAMIP_hist-aer', etc.
% end
% obs_str='SW calc';
% ACSIS_Robson_paper_choose_clims_etc %This sets the ylab, etc.
% 
% %ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic
% 
% 
% 
% 
% %% Nd and CF constant
% sw_str = 'CF and Nd constant';
% clear dat_ukesm_Nd_cf
% dat_ukesm_Nd_cf.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
% dat_ukesm_Nd_cf.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
% dat_ukesm_Nd_cf.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
% dat_ukesm_Nd_cf.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
% dat_ukesm_Nd_cf.years_ukesm_1d = years_sw_calc';
% 
% dat_ukesm_Nd_cf.dat_annual = f_sw_calc * SWTOA_calc_model_annual_Nd_cf;
% %Just repeat with dummy data for ensemble values for now - could actually
% %calculate the individual ensembles.
% a = repmat(f_sw_calc * SWTOA_calc_model_annual_Nd_cf,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm_Nd_cf.dat_annual_ens=a;
% 
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_contributions


% %% LWP constant
% sw_str = 'LWP constant';
% clear dat_ukesm_lwp
% dat_ukesm_lwp.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
% dat_ukesm_lwp.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
% dat_ukesm_lwp.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
% dat_ukesm_lwp.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
% dat_ukesm_lwp.years_ukesm_1d = years_sw_calc';
% 
% dat_ukesm_lwp.dat_annual = f_sw_calc * SWTOA_calc_model_annual_lwp;
% %Just repeat with dummy data for ensemble values for now - could actually
% %calculate the individual ensembles.
% a = repmat(f_sw_calc * SWTOA_calc_model_annual_lwp,[1 1 1 2]);
% a=permute(a,[4 1 2 3]);
% dat_ukesm_lwp.dat_annual_ens=a;
% 
% dat_ukesm_sw_sens = dat_ukesm_lwp;
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_contributions


end

end