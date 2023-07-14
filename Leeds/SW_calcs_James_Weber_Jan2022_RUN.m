% SW_calcs_James_Weber_Jan2022_RUN.m
% Based on ACSIS_Robson_paper_offline_SW_calcs_Sep2020_CALCs.m

%whether to save the resultign aci forcing and where to save it
isave = 1;
save_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/';
save_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/James_King_CESM/forests_UKESM/';

iload = 1;

%% Load the input data from James's NetCDFs.
if iload==1
    SW_calcs_James_Weber_Jan2022_load_input_data
end

%% Run some setup of methods, etc.
SW_calcs_James_Weber_Jan2022_setup

%% Calculate annual means of the true/actual model values

dNd = dat_PD.Nd - dat_PI.Nd;
dNd_3d = dat_PD.Nd_ic - dat_PI.Nd_ic;

%Calculate the annual SWTOA for the SW from model output (i.e., the true
%value we are aiming for).
clear SWTOA_clear_sky_model_annual_PI SWTOA_clear_sky_model_annual_PD SWTOA_model_annual_PI SWTOA_model_annual_PD
clear dNd_annual dNd_3d_annual
for iy=1:size(dat_PI.SW_up_TOA,1)/12 %assuming monthly here
    istart=(iy-1)*12+1;
    SWTOA_clear_sky_model_annual_PI(iy,:,:) = meanNoNan(dat_PI.SW_up_TOA_clear_sky(istart:istart+11,:,:),1);
    SWTOA_clear_sky_model_annual_PD(iy,:,:) = meanNoNan(dat_PD.SW_up_TOA_clear_sky(istart:istart+11,:,:),1);
    SWTOA_model_annual_PI(iy,:,:) = meanNoNan(dat_PI.SW_up_TOA(istart:istart+11,:,:),1);
    SWTOA_model_annual_PD(iy,:,:) = meanNoNan(dat_PD.SW_up_TOA(istart:istart+11,:,:),1);
    dNd_annual(iy,:,:) = meanNoNan(dNd(istart:istart+11,:,:),1);
    dNd_3d_annual(iy,:,:,:) = meanNoNan(dNd_3d(istart:istart+11,:,:,:),1);
end

%% Do the SW calculations for various things
time_mean=0;
time_mean=1;
cf_min = 0.05;

% PI clouds and PI Nd.

cf = dat_PI.totCF;
nd = dat_PI.Nd;
lwp = dat_PI.lwp_ic;
sw_clear = dat_PI.SW_up_TOA_clear_sky;
sw_in = dat_PI.SW_down_TOA;

[SWTOA_calc_model,SWTOA_calc_model_annual,SWTOA_clear_sky_calc_model,SWTOA_clear_sky_calc_model_annual] = SW_calcs_James_Weber_Jan2022_sw_calc(sw_clear,f_upscatter,trans_clear_sky,cf,lwp,nd,sw_in,surf_albedo,tr_for_SW_bias,tr_atmos_constant_val,fac_SWin,cf_min);


%% PD clouds and PD Nd.
cf = dat_PD.totCF;
nd = dat_PD.Nd;
lwp = dat_PD.lwp_ic;
sw_clear = dat_PD.SW_up_TOA_clear_sky;
sw_in = dat_PD.SW_down_TOA;

[SWTOA_calc_model_PD,SWTOA_calc_model_annual_PD,SWTOA_clear_sky_calc_model_PD,SWTOA_clear_sky_calc_model_annual_PD] = SW_calcs_James_Weber_Jan2022_sw_calc(sw_clear,f_upscatter,trans_clear_sky,cf,lwp,nd,sw_in,surf_albedo,tr_for_SW_bias,tr_atmos_constant_val,fac_SWin,cf_min);



            
%% PI clouds with PD Nd.
cf = dat_PI.totCF;
nd = dat_PD.Nd;
lwp = dat_PI.lwp_ic;
sw_clear = dat_PI.SW_up_TOA_clear_sky;
sw_in = dat_PI.SW_down_TOA;

[SWTOA_calc_model_PI_PD_Nd,SWTOA_calc_model_annual_PI_PD_Nd,SWTOA_clear_sky_calc_model_PI_PD_Nd,SWTOA_clear_sky_calc_model_annual_PI_PD_Nd] = SW_calcs_James_Weber_Jan2022_sw_calc(sw_clear,f_upscatter,trans_clear_sky,cf,lwp,nd,sw_in,surf_albedo,tr_for_SW_bias,tr_atmos_constant_val,fac_SWin,cf_min);

%% PD clouds with PI Nd.
cf = dat_PD.totCF;
nd = dat_PI.Nd;
lwp = dat_PD.lwp_ic;
sw_clear = dat_PD.SW_up_TOA_clear_sky;
sw_in = dat_PD.SW_down_TOA;

[SWTOA_calc_model_PD_PI_Nd,SWTOA_calc_model_annual_PD_PI_Nd,SWTOA_clear_sky_calc_model_PD_PI_Nd,SWTOA_clear_sky_calc_model_annual_PD_PI_Nd] = SW_calcs_James_Weber_Jan2022_sw_calc(sw_clear,f_upscatter,trans_clear_sky,cf,lwp,nd,sw_in,surf_albedo,tr_for_SW_bias,tr_atmos_constant_val,fac_SWin,cf_min);


%% Calculate a forcing
aci_PIPD = SWTOA_calc_model_annual_PI_PD_Nd - SWTOA_calc_model_annual;
aci_PDPI = SWTOA_calc_model_annual_PD - SWTOA_calc_model_annual_PD_PI_Nd;
aci = 0.5* (aci_PIPD + aci_PDPI);

%monthly aci too
aci_PIPD_monthly = SWTOA_calc_model_PI_PD_Nd - SWTOA_calc_model;
aci_PDPI_monthly = SWTOA_calc_model_PD - SWTOA_calc_model_PD_PI_Nd;
ERF_Nd_TOA_up_monthly = 0.5* (aci_PIPD_monthly + aci_PDPI_monthly);

%% Save aci forcing
if isave==1
    save_name = [save_dir 'ERF_Nd_' model_type '.mat'];
    lat2d = dat_PI.lat2d;
    lat2d_edges = dat_PI.lat2d_edges;
    lon2d = dat_PI.lon2d;
    lon2d_edges = dat_PI.lon2d_edges;
    ERF_Nd_TOA_up_time_mean = meanNoNan(aci,1);
    Nd_baseline_monthly = dat_PI.Nd;
    Nd_maxforest_monthly = dat_PD.Nd;
    save(save_name,'ERF_Nd_TOA_up_monthly','ERF_Nd_TOA_up_time_mean','lat2d','lat2d_edges','lon2d','lon2d_edges',...
        'Nd_baseline_monthly','Nd_maxforest_monthly','-V7.3');
    mat2nc_Dan(save_name,[save_name '.nc']);
end



%% Maps of various things
i_plot_maps=1;
if i_plot_maps==1    
    SW_calcs_James_Weber_Jan2022_PLOTS        
end





