%can use this to save a specified list of variables witha  tag on the end
%of each
% See case 153 for example of loading and plotting the data

save_or_load = 'save';

tag = 'AMSREvsCF';
tag = 'MODIS_MOD35_15_percent_vsCF';
tag = 'MODIS_MOD06_vsCF';
%tag = 'MODIS_MOD35_vsCF';

tag = 'het_low_SZA';
%tag = 'het_high_SZA';

%tag = 'lowSZA';
%tag = 'highSZA';

%tag = 'Jan'
%tag = 'March';
%tag = 'Sep'

%tag = gcm_str;

%tag='';


filename = ['saved_vars_for_Nd_vs_height_60-70S_60_160W']; 
%filename = ['saved_vars_for_high_minus_low_sza_vs_ctt']; %re
%filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Nd']; %Nd
%filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau']; %tau
filename = ['saved_vars_for_high_minus_low_sza_vs_ctt_Tau_for_ACP']; %tau
%filename = ['saved_vars_for_high_minus_low_sza_vs_gamma_tau']; %re
% filename = ['saved_vars_for_LWP_vs_CF']; 
% filename = ['saved_vars_for_CF_cumulative_PDFs_for_LWP.LT.1']; 
% %filename = ['saved_vars_for_CF_PDFs_for_LWP.LT.1']; 
% filename = ['saved_vars_for_LWP_PDFs_for_CF.GT.0.05']; 
% filename = ['saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1']; 
% filename = ['saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.GT.1']; 
% filename = ['saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.TLWP.LT.1']; 
% %filename = ['saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1']; 
% filename = ['saved_vars_for_TLWP_PDFs_for_LWP.LT.1.AND.PRECIP.GT.0.1']; 
% filename = ['saved_vars_for_LWP_PDFs_for_CF.GT.0.05_all_times_and_locations']; 
% filename = ['saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.LT.1_all_times_and_locations']; 
% %filename = ['saved_vars_for_precip_PDFs_for_CF.GT.0.05.AND.LWP.GT.1_all_times_and_locations']; 
% filename = ['saved_vars_for_fraction_points_95percent_LWP_removed_mphys_all_points_all_times_and_locations']; 
% filename = ['Het_params_comparison_SZA_paper_Cahaan_vs_sigCTT'];
 filename = ['Het_params_comparison_SZA_paper_Cahaan_vs_sigCTT_binned_by_sigCTT'];

%filedir_savevars = '/home/disk/eos1/d.grosvenor/';
%changed to this (from the above) to keep tidier :-
filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
filename_savevars = [filedir_savevars filename '.mat'];

ivar=1; clear vars_to_save


vars_to_save{ivar}='xdat'; ivar=ivar+1;
vars_to_save{ivar}='ydat'; ivar=ivar+1;
vars_to_save{ivar}='errordatU'; ivar=ivar+1;
vars_to_save{ivar}='errordatL'; ivar=ivar+1;
vars_to_save{ivar}='thresh_str'; ivar=ivar+1;
vars_to_save{ivar}='xlab'; ivar=ivar+1;
vars_to_save{ivar}='ylab'; ivar=ivar+1;
vars_to_save{ivar}='titlenam'; ivar=ivar+1;
vars_to_save{ivar}='labs'; ivar=ivar+1;
vars_to_save{ivar}='Ndatap'; ivar=ivar+1;

%run the save script.
save_vars_mat