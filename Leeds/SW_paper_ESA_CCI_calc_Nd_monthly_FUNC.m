function SW_paper_ESA_CCI_monthly_data_read_FUNC()
%Read data from the ESA CCI cloud dataset and save

dat_dir = '/home/disk/eos8/d.grosvenor/ESA_Cloud_CCI/AVHRR_PMv3_L3C_Monthly/';

%var_list={'cer_liq','cot_liq','ctt'};
load_file = [dat_dir 'ESA_Cloud_CCI_Monthly_vars_for_Nd.mat'];
save_file = load_file;
load_file2 = [dat_dir 'ESA_Cloud_CCI_Monthly_Cloud_Fraction.mat'];

%Make the blank arrays
% for ivar=1:length(var_list)
%     eval_str = [var_list{ivar} '=NaN*ones([nt length(lat) length(lon)]);'];
%     eval(eval_str);
% end

dat_cci = load(load_file,'cot_liq','cer_liq','ctt');
dat_cci2 = load(load_file2,'cfc_day','cfc_low');

Wflag='calc';
fad=1;
k=0.8;
[N,H,W,k,Q,cw]=MODIS_N_H_func(dat_cci.cot_liq,dat_cci.cer_liq*1e-6,Wflag,NaN,dat_cci.ctt,fad,k);

%Do some filtering (although more difficult with monthly mean data)

iT268 = find(dat_cci.ctt < 268); %filter for cold/high clouds
N_T268 = N;
N_T268(iT268) = NaN;

icot5 = find(dat_cci.cot_liq<5);
N_T268_cot5 = N_T268;
N_T268_cot5(icot5) = NaN;

ilowCF80 = find(dat_cci2.cfc_low<0.8);
ilowCF20 = find(dat_cci2.cfc_low<0.2);
ilowCF40 = find(dat_cci2.cfc_low<0.4);
N_T268_cf80 = N_T268;
N_T268_cf80(ilowCF80) = NaN;
N_T268_cf80_cot5 = N_T268_cot5;
N_T268_cf80_cot5(ilowCF80) = NaN;

ire3 = find(dat_cci.cer_liq<3);
N_T268_cf80_cot5_re3 = N_T268_cf80_cot5;
N_T268_cf80_cot5_re3(ire3) = NaN;
N_T268_cot5_re3 = N_T268_cot5;
N_T268_cot5_re3(ire3) = NaN;
N_T268_cf20_cot5_re3 = N_T268_cot5_re3;
N_T268_cf20_cot5_re3(ilowCF20) = NaN;
N_T268_cf40_cot5_re3 = N_T268_cot5_re3;
N_T268_cf40_cot5_re3(ilowCF40) = NaN;
save(save_file,'-V7.3','-APPEND','N','H','W','k','Q','cw','N_T268','N_T268_cot5','N_T268_cf80',...
    'N_T268_cf80_cot5','N_T268_cf80_cot5_re3','N_T268_cot5_re3','N_T268_cf20_cot5_re3','N_T268_cf40_cot5_re3');
%,'sw_up_toa','sw_net_coarse','M_coarse_grain','N_coarse_grain','yr_start_UM','yr_end_UM',...
%    'sw_in_TOA_UM','sw_deepc','time_matlab','lat2d_deepc','lon2d_deepc');




