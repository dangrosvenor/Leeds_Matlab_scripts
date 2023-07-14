
dirUM = [UM_base_dir um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon,var_UM,opts);
%clear opts; opts.dummy=NaN;
if ~exist('icheck_time') | icheck_time==1
    %if isfield('opts','time_inds_override
    opts.time_out_check = time_out;
end
clear icheck_time
[time_inds_var,time_out_var] = UM_get_time_inds_from_dat_global(dat_global,time_choice,opts);
eval_str =['[' var_name_out ',nT] = UM_get_time_data_mat_nc(dat_global,time_inds_var,load_type,gcm_Plat2D_UM);']; eval(eval_str);