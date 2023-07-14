file_dir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/';

models = {};
%models{end+1} = {'BCC-ESM1'};
% models{end+1} = {'CNRM-CM6-1'};
% models{end+1} = {'CNRM-ESM2-1'};
 %models{end+1} = {'EC-Earth3-AerChem'};
% %models{end+1} = {'MPI-ESM-1-2-HAM'};
% models{end+1} = {'IPSL-CM6A-LR'};
% models{end+1} = {'MIROC6'};
% models{end+1} = {'MIROC-ES2L'};
% %models{end+1} = {'MRI-ESM2-0'};
% models{end+1} = {'CESM2'};
% models{end+1} = {'CESM2-FV2'};
% models{end+1} = {'CESM2-WACCM'};
% models{end+1} = {'CESM2-WACCM-FV2'};
% models{end+1} = {'NorESM2-MM'};
 models{end+1} = {'NIMS-KMAUKESM1-0-LL'};
% models{end+1} = {'GFDL-CM4'};
% models{end+1} = {'GFDL-ESM4'};


var = 'Nd_clw_weighted_ESGF_no_dz';

for i=1:length(models)
   cmip6_str = models{i}{1};
   mat_file = [file_dir 'EGSF_ensemble_timeseries_' cmip6_str '_all_' var '.mat'];
   mat_obj = matfile(mat_file);
   qpcolor(meanNoNan(mat_obj.dat_annual,1)/1e6);
   title(cmip6_str);
   caxis([0 500]);
end

%%
models = {'HADGEM3_GC31_LL'};

var = 'Nd_clw_weighted_ESGF';

for i=1:length(models)
   cmip6_str = models{i};
   mat_file = [file_dir 'EGSF_ensemble_timeseries_' cmip6_str '_all_' var '.mat'];
   mat_obj = matfile(mat_file);
   qpcolor(meanNoNan(mat_obj.dat_annual,1)/1e6);
   title(cmip6_str);
   caxis([0 500]);
end