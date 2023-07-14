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
var_UM_list = {'rsut','clt'};
%var_UM_list = {'rsut','clt','od550dust','od865dust'};
%var_UM_list = {'psl','pr'}; %prra and tos are from the ocean model at higher res with different lat lon names.
 %Doesn't have the correct merged files for psl
%var_UM_list = {'psl'}; %prra and tos are from the ocean model at higher res with different lat lon names.
 %Doesn't have the correct merged files for psl
%var_UM_list = {'scldncl'}; %prra and tos are from the ocean model at higher res with different lat lon names.
 %Doesn't have the correct merged files for psl
 
 
%var_UM_list = {'ccldncl','scldncl','ps'}; %The first two here (convective and stratiform top CDNC) have data every day and it causes
%it to run out of memory. Need to process into monthly first - or make the
%Nd fields myself on Jamsin from the monthly data. 
% ps is not present - need to copy across from Jasmin.

%N.B. - the cloud masks for calipso cosp variables it automatically
%dealt with in the ACSIS_dat_trends_load_ensemble function

for ivar=1:length(var_UM_list)
    
    var_UM = var_UM_list{ivar};
    
    ACSIS_dat_trends_load_ensemble_esgf_HadGEM(var_UM);

end

