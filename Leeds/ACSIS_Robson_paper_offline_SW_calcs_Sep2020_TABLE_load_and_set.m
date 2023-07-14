
%var_str='Albedo';
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_SW_TOA_up_calc_for_region_4__ocean_only,_no_sea-ice_' model_str 'albedo_vary__' period_str '.mat'];

switch model_str
    case 'AerChemMIP_Aerosol_Proxy_'
        tr_file = remove_character(tr_file,model_str,'UKESM1-AerChemMIP_control_');
end

sw_trends = load(tr_file);
%Full SW trend values from the model
eval_str=['table_vals.' var_str period_str_char ' = fdT_' period_str '*fscale*sw_trends.trend_dat_box{1,itr}.coeffs(2);']; eval(eval_str);
eval_str=['table_vals.' var_str period_str_char 'un = fdT_' period_str '*fscale*sw_trends.trend_dat_box{1,itr}.uncer_max;']; eval(eval_str);

 
switch model_str
    case 'AerChemMIP_Aerosol_Proxy_'
        tr_file = remove_character(tr_file,'UKESM1-AerChemMIP_control_','AerChemMIP_hist-piAer_');
        
        sw_trends = load(tr_file);
        %Do the full AerChemMIP minus the piAer to estimate the aerosol
        %effect
        eval_str=['var = fdT_' period_str '*fscale*sw_trends.trend_dat_box{1,itr}.coeffs(2);']; eval(eval_str);
        eval_str=['var_un = fdT_' period_str '*fscale*sw_trends.trend_dat_box{1,itr}.uncer_max;']; eval(eval_str);
        
        eval_str=['table_vals.' var_str period_str_char ' = table_vals.' var_str period_str_char ' - var;']; eval(eval_str);
        eval_str=['table_vals.' var_str period_str_char 'un = table_vals.' var_str period_str_char 'un + var_un;']; eval(eval_str);

end

