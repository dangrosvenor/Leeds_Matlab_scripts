function [LES_struct] = calc_mean_fields_Zhibo_FUNCS_timemean(LES_struct_in,var_str_in,var_str2_in,rt_str,itime,ntimes)

LES_struct = LES_struct_in;
var_str = [var_str2_in '_' var_str_in];
var_str2 = var_str2_in;

if itime==1
    LES_struct.(['timemean_' var_str2 '_' rt_str '_vs_res']) = 0;
end
LES_struct.(['timemean_' var_str2 '_' rt_str '_vs_res']) =  LES_struct.(['timemean_' var_str2 '_' rt_str '_vs_res']) + LES_struct.(['mean_' var_str2 '_' rt_str '_vs_res']) ./ ntimes; %save running average for time mean


