expt_str = 'AerChemMIP_hist-piAer'; expt_str2='AerChemMIP_piaer';

switch var_ukesm
    case 'scldncl'
        var_DAMIP = 'Nd_clw_weighted_ESGF';
        fscale_DAMIP = 1e-6; %convert from m^-3 to cm^-3
    case 'Nd_clw_weighted_ESGF'
        var_DAMIP = 'Nd_clw_weighted_ESGF_total_column_to_zdomain_top';
        fscale_DAMIP = dat_ukesm.fscale;
    otherwise
        var_DAMIP = var_ukesm;
        fscale_DAMIP = dat_ukesm.fscale;
end

load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str '_all_' var_DAMIP '.mat'];
load_file_PI = load_file;

eval_str =['[dat_PI_' expt_str2 ',dat_ukesm_' expt_str2 ',dat_ukesm_' expt_str2 '_DJF,dat_ukesm_' expt_str2 ...
    '_JJA,dat_ukesm_' expt_str2 '_MAM,dat_ukesm_' expt_str2 '_SON]='...
    'ACSIS_Robson_paper_load_data_generic_FUNC(load_file,load_file_PI,var_DAMIP,fscale_DAMIP);'];
eval(eval_str);

%Run script to choose the data for the selected season and then to average
%over the box region.
yr_start_trend_box2 = yr_start_trend_box; yr_end_trend_box2 = yr_end_trend_box;
end_str = ['_' expt_str2];
ACSIS_Robson_paper_box_means

%%
