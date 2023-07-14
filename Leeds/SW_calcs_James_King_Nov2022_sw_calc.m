function [SWTOA_calc_model,SWTOA_calc_model_annual,SWTOA_clear_sky_calc_model,SWTOA_clear_sky_calc_model_annual] = SW_calcs_James_King_Nov2022_sw_calc(SW_up_TOA_clear_sky,f_upscatter,trans_clear_sky,cf,lwp,nd,sw,surf_albedo,tr_for_SW_bias,tr_atmos_constant_val,fac_SWin,cf_min,i_Liu);
%nd in cm^{-3}, W in kg/m3, cf between 0 and 1.


%cf = lowcf_dat_ens_mean;
%cf = totcf_dat_ens_mean;

% cf = dat_PI.tot_CF;
% nd = dat_PI.Nd;
% lwp = dat_PI.lwp_ic;

switch tr_for_SW_bias
    case 'constant value'
        trans_orig = tr_atmos_constant_val; %if constant value
        %sw = meanNoNan(SW_in,3);
        % sw = fac_SWin * SW_down_TOA_dat_ens_mean .* trans_orig;
        sw = fac_SWin * sw;
        
    case 'time mean from model'
        trans_orig = tr_atmos_mean;
        %sw = meanNoNan(SW_in,3);
        sw = fac_SWin * SW_down_TOA_dat_ens_mean .* trans_orig;
        
    case 'time varying monthly'
        %f_trans = sqrt(1.21);
        f_trans = 1;
        %f_trans = 0.95;
        f_trans = 1.25;
        trans_orig = f_trans * transmission_atmos_down; %can only really use the downwards transmission since upwards includes
        %the scattered light from the sun and not just that originating
        %from the surface.
        
        trans_orig(trans_orig>1) = 1;
        
        %sw = fac_SWin * sw_surf_down_clear_dat_ens_mean; %SW_down_TOA_dat_ens_mean .* transmission_atmos_down;
        %sw = fac_SWin * SW_down_TOA_dat_ens_mean .* trans_orig;
        sw = fac_SWin * SW_down_TOA_dat_ens_mean; %Now running calc_SW2 that uses TOA SW.
        
    case 'time varying annual mean'
        
        clear trans_orig_ann
        for iy=1:length(years_sw_calc)
            istart=(iy-1)*12+1;
            %trans_orig_ann(iy,:,:) = meanNoNan(transmission_atmos_down(istart:istart+11,:,:),1,'mean',1,1,SW_down_TOA_dat_ens_mean(istart:istart+11,:,:));
            trans_orig_ann(iy,:,:) = 1.27*max(transmission_atmos_down(istart:istart+11,:,:),[],1);
        end
        trans_orig = NaN*ones(size(transmission_atmos_down));
        for im=1:size(transmission_atmos_down,1)
            iy = ceil(im/12);
            trans_orig(im,:,:) = trans_orig_ann(iy,:,:);
        end
        
        %sw = fac_SWin * sw_surf_down_clear_dat_ens_mean; %SW_down_TOA_dat_ens_mean .* transmission_atmos_down;
        sw = fac_SWin * SW_down_TOA_dat_ens_mean .* trans_orig;
end


str_type_calc = 'time-mean';


clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=i_Liu; opts.iprovide_SWclear=1; opts.SWclear = SW_up_TOA_clear_sky;
%[Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model,SWTOA_clear_sky_calc_model] = calc_SW2(cf,lwp*1e3,nd,sw,surf_albedo,trans_orig,cf_min,0,NaN,i_Liu);
[Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model,SWTOA_clear_sky_calc_model] = calc_SW3(opts,NaN,NaN,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,surf_albedo,trans_orig);

SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);

clear SWTOA_calc_model_annual SWTOA_clear_sky_calc_model_annual
for iy=1:size(SWTOA_calc_model,1)/12 %assuming monthly here
    istart=(iy-1)*12+1;
    SWTOA_calc_model_annual(iy,:,:) = meanNoNan(SWTOA_calc_model(istart:istart+11,:,:),1);
    SWTOA_clear_sky_calc_model_annual(iy,:,:) = meanNoNan(SWTOA_clear_sky_calc_model(istart:istart+11,:,:),1);
end

% clear surf_albedo_annual
% for iy=1:length(years_sw_calc)
%     istart=(iy-1)*12+1;
%     surf_albedo_annual(iy,:,:) = meanNoNan(surf_albedo(istart:istart+11,:,:),1);
% end



