
imodel=1;
clear MIP_multi

%MIP_multi{imodel} = 'CMIP6'; imodel=imodel+1; %UKESM1
%MIP_multi{imodel} = 'AMIP'; imodel=imodel+1;
%MIP_multi{imodel} = 'PI'; imodel=imodel+1;
%MIP_multi{imodel} = 'DAMIP'; imodel=imodel+1;
%MIP_multi{imodel} = 'HADGEM3_GC31_LL'; imodel=imodel+1;
%MIP_multi{imodel} = 'Nudged'; imodel=imodel+1;
MIP_multi{imodel} = 'AerChemMIP'; imodel=imodel+1;

 %MIP_multi{imodel} = 'CMIP6 HadGEM3-GC31-LL'; imodel=imodel+1;
 %MIP_multi{imodel} = 'CMIP6 UKESM'; imodel=imodel+1;
%MIP_multi{imodel} = 'CMIP6 BCC-ESM1'; imodel=imodel+1;
%MIP_multi{imodel} = 'CMIP6 CNRM-CM6-1'; imodel=imodel+1;
%MIP_multi{imodel} = 'CMIP6 CNRM-ESM2-1'; imodel=imodel+1;
% % MIP_multi{imodel} = 'CMIP6 EC-Earth3-AerChem'; imodel=imodel+1;
%  MIP_multi{imodel} = 'CMIP6 MPI-ESM-1-2-HAM'; imodel=imodel+1; %Only 60 timesteps
%   MIP_multi{imodel} = 'CMIP6 IPSL'; imodel=imodel+1;
   %MIP_multi{imodel} = 'CMIP6 MIROC6'; imodel=imodel+1;
%   MIP_multi{imodel} = 'CMIP6 MIROC-ES2L'; imodel=imodel+1;
  %MIP_multi{imodel} = 'CMIP6 MRI-ESM2-0'; imodel=imodel+1; %Only 60 timesteps
%  MIP_multi{imodel} = 'CMIP6 CESM2'; imodel=imodel+1; %reff
  %MIP_multi{imodel} = 'CMIP6 CESM2-FV2'; imodel=imodel+1;
  %MIP_multi{imodel} = 'CMIP6 CESM2-WACCM'; imodel=imodel+1;
  %MIP_multi{imodel} = 'CMIP6 CESM2-WACCM-FV2'; imodel=imodel+1;
% % % MIP_multi{imodel} = 'CMIP6 NorESM2-MM';imodel=imodel+1;
   %MIP_multi{imodel} = 'CMIP6 NIMS-KMA UKESM1-0-LL'; imodel=imodel+1;
% MIP_multi{imodel} = 'CMIP6 GFDL-CM4'; imodel=imodel+1;
%    MIP_multi{imodel} = 'CMIP6 GFDL-ESM4'; imodel=imodel+1;
%  MIP_multi{imodel} = 'CMIP6 INM-CM4-8'; imodel=imodel+1;
%  MIP_multi{imodel} = 'CMIP6 INM-CM5-0'; imodel=imodel+1;

%expt is currently used for DAMIP, nudged and AerChemMIP only
%DAMIP
%expt='hist-aer';
%expt='hist-GHG';
%expt='hist-nat';

expt = 'hist-piAer';

%Nudged run experiment
%expt='u-by844'; run_str='nudged_WIND_only';
%expt='u-bz785'; run_str='nudged_WIND_only_abv_BL';

for im=1:length(MIP_multi)
    MIP = MIP_multi{im};
    ACSIS_dat_trends_load_ESGF_ensemble_multi_vars    
end


