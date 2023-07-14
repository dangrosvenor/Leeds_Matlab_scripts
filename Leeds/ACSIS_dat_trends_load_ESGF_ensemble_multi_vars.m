%This can be run from ACSIS_load_ESGF_ensemble_multi_vars_MULTI_model.m to
%run for multiple models.

%function []=ACSIS_dat_trends_load_ensemble_multi_vars()
% Loads in the UM .nc files for the ensembles, does some processing and
% stores as .mat files.

expt='';

MIP = 'CMIP6'; %some of the UKESM1 variables are set to use this and others CMIP6 UKESM1 - it depends on where they are stored
%Check in the ACSIS_dat_trends_load_ensemble_esgf script.
%MIP = 'CMIP6 UKESM'; %UKESM1
%MIP = 'AMIP';
%MIP = 'PI';
MIP = 'DAMIP'; 
MIP = 'HADGEM3_GC31_LL';
%MIP = 'Nudged';
%MIP = 'AerChemMIP';
MIP = 'CMIP6 UKESM1-AerChemMIP_control';

%MIP = 'CMIP6 BCC-ESM1';
 %MIP = 'CMIP6 CNRM-CM6-1';
 %MIP = 'CMIP6 CNRM-ESM2-1';
% MIP = 'CMIP6 EC-Earth3-AerChem';
 %MIP = 'CMIP6 MPI-ESM-1-2-HAM'; %Only 60 timesteps
 %MIP = 'CMIP6 IPSL';
 %MIP = 'CMIP6 MIROC6';
% MIP = 'CMIP6 MIROC-ES2L';
 %MIP = 'CMIP6 MRI-ESM2-0'; %Only 60 timesteps
 %MIP = 'CMIP6 CESM2'; %
% MIP = 'CMIP6 CESM2-FV2';
% MIP = 'CMIP6 CESM2-WACCM';
% MIP = 'CMIP6 CESM2-WACCM-FV2';
% MIP = 'CMIP6 NorESM2-MM';
% MIP = 'CMIP6 NIMS-KMA UKESM1-0-LL';
% MIP = 'CMIP6 GFDL-CM4';
 %MIP = 'CMIP6 GFDL-ESM4';
 %MIP = 'CMIP6 INM-CM4-8';
 %MIP = 'CMIP6 INM-CM5-0';

%expt is currently used for DAMIP, nudged and AerChemMIP only
%DAMIP
 expt='hist-aer';
 expt='hist-GHG';
 expt='hist-nat';

%expt = 'hist-piAer'; %For AerChemMIP

%Nudged run experiment
%expt='u-by844'; run_str='nudged_WIND_only';
%expt='u-bz785'; run_str='nudged_WIND_only_abv_BL';



%var_UM = 'Nd_cf_weighted_UKESM_ztop'; 
%var_UM = 'Nd_cf_weighted_UKESM';

%var_UM_list = {'calipso_total_cloud_amount','SW_up_TOA','calipso_low_cloud_amount','Nd_cf_weighted_UKESM','LWP_2-391','calipso_mid_cloud_amount'};
    %,'calipso_high_cloud_amount'}; %data missing for high at moment
%var_UM_list = {'Nd_cf_weighted_UKESM'};
%var_UM_list = {'calipso_total_cloud_amount'};
%var_UM_list = {'SW_down_TOA'};
%var_UM_list = {'clt','rsut','ts'};
%var_UM_list = {'ts'};
%var_UM_list = {'rsut','clt','od550dust','od865dust'};
% var_UM_list = {'psl','pr'}; %prra and tos are from the ocean model at higher res with different lat lon names.
%  %Doesn't have the correct merged files for psl
% var_UM_list = {'psl'}; %prra and tos are from the ocean model at higher res with different lat lon names.
%  %Doesn't have the correct merged files for psl
% var_UM_list = {'scldncl'}; %This is daily data, so need to process into monthly data first using ACSIS_Robson_paper_calc_monthly_Nd_ESGF.m
% %and then run this script
% var_UM_list = {'lwp','lwpic'};
% %var_UM_list = {'rsdt'}; %downwelling TOA SW
% %var_UM_list = {'clisccp'};
% var_UM_list = {'lwp'};
% var_UM_list = {'LWP_2-391'};
% var_UM_list = {'IWP_2-392','tot_cloud_amount_in_rad'};
% %var_UM_list = {'tot_cloud_amount_in_rad'};
% %var_UM_list = {'lwpic'}; %in-cloud LWP (divided by the monthly clt).
% %var_UM_list = {'clt'};
% var_UM_list = {'rsut'};
%var_UM_list = {'rsut','clt'};
%var_UM_list = {'SW_up_TOA'};
%var_UM_list = {'Nd_clw_weighted_ESGF'};
%var_UM_list = {'Nd_clw_weighted_ESGF_no_dz'};
%var_UM_list = {'Nd_clw_weighted_ESGF_no_dz_div_CF_total_column_to_zdomain_top'}; %
%var_UM_list = {'Nd_clw_weighted_ESGF_no_dz_div_CF'};
%var_UM_list = {'Nd_clw_weighted_ESGF_no_dz_div_CF_no_ice_total_column_to_zdomain_top'}; %
%var_UM_list = {'Nd_clw_weighted_ESGF_no_dz_no_ice_total_column_to_zdomain_top'}; %

%var_UM_list = {'rsds'}; %Surface downwelling SW.
%var_UM_list = {'prw'}; %Surface downwelling SW.
%var_UM_list = {'rsutcs'};
%var_UM_list = {'rsuscs','rsdscs'};
var_UM_list = {'od550aer'};
%var_UM_list = {'od550aerso'}; %Not available for DAMIP - only HADGEM and UKESM1
%var_UM_list = {'abs550aer'};
var_UM_list = {'od550dust'};
var_UM_list = {'od550tot'}; %od550aer + od550dust
%var_UM_list = {'od550aer','abs550aer'};
%var_UM_list = {'reffclwtop'};
%var_UM_list = {'cldnvi'};

%var_UM_list = {'Nd_clw_weighted_ESGF_total_column_to_zdomain_top','clt','lwpic','rsutcs','rsuscs','rsdscs','rsdt'};
%var_UM_list = {'Nd_clw_weighted_ESGF','clt','lwpic','rsutcs','rsuscs','rsdscs','rsdt'};
%var_UM_list = {'rsuscs','rsdscs'};

%N.B. - the cloud masks for calipso cosp variables it automatically
%dealt with in the ACSIS_dat_trends_load_ensemble function

for ivar=1:length(var_UM_list)
    
    var_UM = var_UM_list{ivar}
    
    ACSIS_dat_trends_load_ensemble_esgf(var_UM,MIP,expt);

end

