% Based on SW_calcs_James_Weber_Jan2022_RUN.m
% which was in turn based on ACSIS_Robson_paper_offline_SW_calcs_Sep2020_CALCs.m

%whether to save the resultign aci forcing and where to save it
isave = 0;
save_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/';



%% Run some setup of methods, etc.
SW_calcs_Hawaii_setup_vals
SW_calcs_Hawaii_setup

%% Sort out the cloud variables
%CF proxy using LWP threshold - either zero or one depending if LWP is
%below or above thereshold.

CF_type = 'LWP_threshold';
CF_type = 'CF_subgrid_vert_int';

switch CF_type
    case 'LWP_threshold'
        
        thresh_lwp_CF = 0.1; %g/m2
        thresh_lwp_CF = 5; %g/m2
        
        %f0
        f0_orig = zeros(size(LWP_PI_ALL)); %f0_orig(inan) = NaN;
        f0_orig(isnan(LWP_PI_ALL)) = NaN;
        icloud = find(LWP_PI_ALL>thresh_lwp_CF);
        f0_orig(icloud) = 1;
        
        W0_orig = NaN*ones(size(LWP_PI_ALL));
        W0_orig(icloud) = LWP_PI_ALL(icloud);
        N0_orig = NaN*ones(size(LWP_PD_ALL));
        N0_orig(icloud) = Nd_PI_ALL(icloud)/1e6; %in per cm3
        
        %f1
        f1_orig = zeros(size(LWP_PD_ALL)); %f1_orig(inan) = NaN;
        f1_orig(isnan(LWP_PD_ALL)) = NaN;
        icloud = find(LWP_PD_ALL>thresh_lwp_CF);
        f1_orig(icloud) = 1;
        
        W1_orig = NaN*ones(size(LWP_PD_ALL));
        W1_orig(icloud) = LWP_PD_ALL(icloud);
        N1_orig = NaN*ones(size(LWP_PD_ALL));
        N1_orig(icloud) = Nd_PD_ALL(icloud)/1e6; %in per cm3


    case 'CF_subgrid_vert_int'
        Hawaii_calc_LWPic_using_subgrid_CF
        
        
        
        
end


[f0,N0,W0,f1,N1,W1,W0_2,N0_2,W1_2,N1_2] = SWTOA_offline_make_cloud_variables(f0_orig,N0_orig,W0_orig,f1_orig,N1_orig,W1_orig,cf_min);



%% Do the SW calculations for various things

%clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1;
clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=0; opts.iprovide_SWclear=1;
%opts.SWclear = SW_up_clean_clear_TOA_PI_ALL; %Probably should use the clear-sky fluxes here


%[Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model,SWTOA_clear_sky_calc_model] = ...
%calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,surf_albedo,tr_atmos_constant_val);

tr_atmos_constant_val = 0.89;
tr_atmos_constant_val = 0.88;
%tr_atmos_constant_val = 0.93;
%tr_atmos_constant_val = 1.0;

SW_in_TOA = SW_down_clean_TOA_PI_ALL;

opts.SWclear = SW_up_clear_TOA_PD_ALL;

%PI baseline calc
%opts.SWclear = SW_up_clear_TOA_PI_ALL;
[Ac,tau,A_calc,SW_PI_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f0,W0,N0,SW_in_TOA,surf_albedo,tr_atmos_constant_val);
[Ac,tau,A_calc,SW_cf_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f1,W0_2,N0_2,SW_in_TOA,surf_albedo,tr_atmos_constant_val);
[Ac,tau,A_calc,SW_lwp_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f0,W1,N0,SW_in_TOA,surf_albedo,tr_atmos_constant_val);
[Ac,tau,A_calc,SW_Nd_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f0,W0,N1,SW_in_TOA,surf_albedo,tr_atmos_constant_val);

%PD baseline calc
%opts.SWclear = SW_up_clear_TOA_PD_ALL; %Should prob use PI clear-sky since
%are just trying to calculated the ACI effects here not direct effects.
[Ac,tau,A_calc,SW_PD_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f1,W1,N1,SW_in_TOA,surf_albedo,tr_atmos_constant_val);
[Ac,tau,A_calc,SW_cf2_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f0,W1_2,N1_2,SW_in_TOA,surf_albedo,tr_atmos_constant_val);
[Ac,tau,A_calc,SW_lwp2_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f1,W0,N1,SW_in_TOA,surf_albedo,tr_atmos_constant_val);
[Ac,tau,A_calc,SW_Nd2_TOA] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,f1,W1,N0,SW_in_TOA,surf_albedo,tr_atmos_constant_val);

%% Estimated forcings

SW_estimated_forcing_tot_TOA = SW_PD_TOA - SW_PI_TOA;
SW_estimated_forcing_tot_TOA_timemean2 = meanNoNan(SW_PD_TOA,3) - meanNoNan(SW_PI_TOA,3);

%Using PI as a baseline
SW_estimated_forcing_cf_TOA = SW_cf_TOA - SW_PI_TOA;
SW_estimated_forcing_Nd_TOA = SW_Nd_TOA - SW_PI_TOA;
SW_estimated_forcing_lwp_TOA = SW_lwp_TOA - SW_PI_TOA;
SW_estimated_forcing_sum_TOA_ALL = SW_estimated_forcing_cf_TOA + SW_estimated_forcing_Nd_TOA + SW_estimated_forcing_lwp_TOA;

%Using PD as a baseline
SW_estimated_forcing_cf_PD_TOA = SW_PD_TOA - SW_cf2_TOA;
SW_estimated_forcing_Nd_PD_TOA = SW_PD_TOA - SW_Nd2_TOA; 
SW_estimated_forcing_lwp_PD_TOA = SW_PD_TOA - SW_lwp2_TOA;
SW_estimated_forcing_sum_PD_TOA_ALL = SW_estimated_forcing_cf_PD_TOA + SW_estimated_forcing_Nd_PD_TOA + SW_estimated_forcing_lwp_PD_TOA;

%Take avearges of PI and PD baseline values and time averages
SW_estimated_forcing_cf_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_cf_PD_TOA,3) + meanNoNan(SW_estimated_forcing_cf_TOA,3) )./2;
SW_estimated_forcing_Nd_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_Nd_PD_TOA,3) + meanNoNan(SW_estimated_forcing_Nd_TOA,3) )./2;
SW_estimated_forcing_lwp_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_lwp_PD_TOA,3) + meanNoNan(SW_estimated_forcing_lwp_TOA,3) )./2;
SW_estimated_forcing_sum_PIPD_TOA_timemean = SW_estimated_forcing_cf_PIPD_TOA_timemean + SW_estimated_forcing_Nd_PIPD_TOA_timemean + SW_estimated_forcing_lwp_PIPD_TOA_timemean;
SW_estimated_forcing_tot_PIPD_TOA_timemean = meanNoNan(SW_estimated_forcing_tot_TOA,3);

%As above, but not using time averages.
SW_estimated_forcing_cf_PIPD_TOA_ALL = ( SW_estimated_forcing_cf_PD_TOA + SW_estimated_forcing_cf_TOA )./2;
SW_estimated_forcing_Nd_PIPD_TOA_ALL = ( SW_estimated_forcing_Nd_PD_TOA + SW_estimated_forcing_Nd_TOA )./2;
SW_estimated_forcing_lwp_PIPD_TOA_ALL = ( SW_estimated_forcing_lwp_PD_TOA + SW_estimated_forcing_lwp_TOA )./2;
SW_estimated_forcing_sum_PIPD_TOA_ALL = SW_estimated_forcing_cf_PIPD_TOA_ALL + SW_estimated_forcing_Nd_PIPD_TOA_ALL + SW_estimated_forcing_lwp_PIPD_TOA_ALL;

ntimes = size(SW_estimated_forcing_cf_PD_TOA,3);
SW_estimated_forcing_cf_PIPD_TOA_timemean_rep = repmat(SW_estimated_forcing_cf_PIPD_TOA_timemean,[1 1 ntimes]);
SW_estimated_forcing_Nd_PIPD_TOA_timemean_rep = repmat(SW_estimated_forcing_Nd_PIPD_TOA_timemean,[1 1 ntimes]);
SW_estimated_forcing_lwp_PIPD_TOA_timemean_rep = repmat(SW_estimated_forcing_lwp_PIPD_TOA_timemean,[1 1 ntimes]);
SW_estimated_forcing_sum_PIPD_TOA_timemean_rep = repmat(SW_estimated_forcing_sum_PIPD_TOA_timemean,[1 1 ntimes]);
SW_estimated_forcing_tot_PIPD_TOA_timemean_rep = repmat(SW_estimated_forcing_tot_PIPD_TOA_timemean,[1 1 ntimes]);
SW_estimated_forcing_tot_PIPD_TOA_timemean2_rep = repmat(SW_estimated_forcing_tot_TOA_timemean2,[1 1 ntimes]);

%% Save aci forcing
if isave==1
    save_name = [save_dir 'ERF_Nd_' model_type '.mat'];
    lat2d = dat_PI.lat2d;
    lat2d_edges = dat_PI.lat2d_edges;
    lon2d = dat_PI.lon2d;
    lon2d_edges = dat_PI.lon2d_edges;
    ERF_Nd_TOA_up = meanNoNan(aci,1);
    save(save_name,'ERF_Nd_TOA_up_monthly','ERF_Nd_TOA_up','lat2d','lat2d_edges','lon2d','lon2d_edges','-V7.3');
    mat2nc_Dan(save_name,[save_name '.nc']);
end



%% Maps of various things
i_plot_maps=0;
if i_plot_maps==1    
    SW_calcs_Hawaii_PLOTS        
end





