%function []=ACSIS_dat_trends_load_ensemble_multi_vars()
% Loads in the UM .nc files for the ensembles, does some processing and
% stores as .mat files.

%var_UM = 'Nd_cf_weighted_UKESM_ztop'; 
%var_UM = 'Nd_cf_weighted_UKESM';

%var_UM_list = {'calipso_total_cloud_amount','SW_up_TOA','calipso_low_cloud_amount','Nd_cf_weighted_UKESM','LWP_2-391','calipso_mid_cloud_amount'};
    %,'calipso_high_cloud_amount'}; %data missing for high at moment
%var_UM_list = {'Nd_cf_weighted_UKESM'};
%var_UM_list = {'calipso_total_cloud_amount'};
%var_UM_list = {'SW_down_TOA'};
%var_UM_list = {'clt','rsut','ts'};
%var_UM_list = {'ts'};
var_UM_list = {'scldncl'}; %Cloud top CDNC - this is daily data, so need to process into monthly data first using ACSIS_Robson_paper_calc_monthly_Nd_ESGF.m
%and then run this script
var_UM_list = {'SW_up_TOA'};
%var_UM_list = {'Nd_cloud_top','LWP_2-391', 'surface_temp', 'tot_cloud_amount_in_rad'};


%N.B. - the cloud masks for calipso cosp variables it automatically
%dealt with in the ACSIS_dat_trends_load_ensemble function

for ivar=1:length(var_UM_list)
    
    var_UM = var_UM_list{ivar};
    
    ACSIS_dat_trends_load_nudged_ACSIS_run(var_UM);

end

