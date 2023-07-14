
switch expt_str
    case 'AerChemMIP_hist-AerProxy'
        expt_str_control = 'UKESM1-AerChemMIP_control';
        load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str_control '_all_' var_DAMIP '.mat'];
        
        expt_str_piAer = 'AerChemMIP_hist-piAer';
        load_file2 = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str_piAer '_all_' var_DAMIP '.mat'];
        
        
    otherwise
        load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str '_all_' var_DAMIP '.mat'];
        load_file2 = '';
end