% SW_calcs_James_Weber_Jan2022_setup
% Called from SW_calcs_James_Weber_Jan2022_RUN.m

%% 2) Set up stuff for the SW TOA calcs from monthly means - likely inaccuracies due to using monthly mean SWTOA down, etc.




trans_clear_sky = 0.95;
trans_clear_sky = 1.0;

%N.B. - can use the clear-sky downwelling at the surface to calc this
%instead of calculating f_upscatter
%f_upscatter = 0.4;

f_up_method = 'sza detailed';
f_up_method = 'use model clear-sky';

switch f_up_method
    case 'use model clear-sky'
        f_upscatter = NaN; %Not needed if using the clear-sky SW
        
    case 'lat'
        
        c = 0.5;
        m = (0.2 - c)/(1-0);
        
        c = 0.8;
        m = (0.4 - c)/(1-0);
        
        f_up= m.*cos(gcm_Plat2D_UM.*pi/180) + c;
        f_upscatter = repmat(f_up,[1 1 size(SW_up_TOA_dat_ens_mean,1)]);
        f_upscatter = permute(f_upscatter,[3 1 2]);
        
        
        
    case 'sza'
        c = 0.6;
        m = (0.2 - c)/(1-0);
        time_LT_hr = 10; time_LT_min = 30; %10:30 LT
        
        c = 0.5;
        m = (0.2 - c)/(1-0);
        time_LT_hr = 9; time_LT_min = 30; %10:30 LT
        
        %Calculate SZA just using the 15th of each month for one year
        %at 0 deg lon.
        
        %mat_times_monthly = datenum(years_sw_calc_monthly,months_sw_calc_monthly,15,time_LT_hr,time_LT_min,0);
        mat_times_monthly = datenum(1850*ones([1 12]),[1:12],15,time_LT_hr,time_LT_min,0);
        mat_times_ALL = repmat(mat_times_monthly',[1 size(SW_up_TOA_dat_ens_mean,2)]);
        lats_ALL = repmat(gcm_Plat2D_UM(:,1),[1 size(mat_times_ALL,1)]); lats_ALL = permute(lats_ALL,[2 1]);
        lons_ALL = zeros(size(mat_times_ALL)); %Use UTC
        sza = sun_pos(mat_times_ALL,lats_ALL,lons_ALL); sza(sza>90) = 90;
        sza_ALL = repmat(sza,[size(SW_up_TOA_dat_ens_mean,1)/12 1 size(SW_up_TOA_dat_ens_mean,3)]);
        
        f_upscatter = m.*cos(sza_ALL.*pi/180) + c;
        
    case 'sza detailed'
        c = 0.6;
        m = (0.2 - c)/(1-0);
        %time_LT_hr = 10; time_LT_min = 30; %10:30 LT
        
        %c = 0.5; m = (0.2 - c)/(1-0);
        c = 0.7; m = (0.0 - c)/(1-0);
        
        
        c = 0.7; m = (0.15 - c)/(1-0.4);
        % c = 0.8; m = (0.05 - c)/(1-0.4);             
        
        
        %Load the SZA data calculated in ACSIS_Robson_paper_offline_SW_calcs_Sep2020_CALCs.m
        dt=0.5;
        loaddir = '/home/disk/eos15/d.grosvenor/UM/ACSIS_Nd_trends/';
        loadname_sza = [loaddir 'sza_' num2str(dt) '_SW_calcs.mat'];
        
        load(loadname_sza,'me_cos'); %[12 144] array of mean cos(sza) for 15th of each month for each UM lat (N96)
            
        %replicate to all years and lons
        siz = size(dat_PI.SW_up_TOA);
        %me_cos_ALL = repmat(me_cos,[size(SW_up_TOA_dat_ens_mean,1)/12 1 size(SW_up_TOA_dat_ens_mean,3)]);
        me_cos_ALL = repmat(me_cos,[siz(1)/12 1 siz(3)]);
            
            
            %save(savename_sza,'sza','me_cos','me_cos_ALL','dt','-V7.3');                   
        
        f_upscatter = m.*me_cos_ALL + c;
                        
end



%Whether to use a constant albedo, or one calculated using surface
%clear-sky fluxes.
albedo_type = 'constant';
albedo_type = 'calculated monthly';
%albedo_type = 'calculated annual mean';


tr_for_SW_bias = 'constant value';
%tr_for_SW_bias = 'time mean from model';
%tr_for_SW_bias = 'time varying monthly';
%tr_for_SW_bias = 'time varying annual mean';


tr_atmos_constant_val = 0.6; %previously used when needed 1.7x factor - use this for now to save re-doing for UKESM
%tr_atmos_constant_val = 0.75;
%tr_atmos_constant_val = 0.83; %Worked better in tests with GC20 ACSIS nudged instantaneous runs.
%tr_atmos_constant_val = 0.7821; %Best match value for DAMIP-Hist-Aer - actually changing trans and f_sw_calc number are equivalent.
%I think. Except that sw_calc is propto trans.^2. So if wanted to
%adjust trans would need to make it. sw1=f*sw0 =
%(trans1/trans0).^2 *sw0  = (trans1/0.83).^2 *sw0. trans1 =
%sqrt(f)*trans0 = sqrt(0.8879)*0.83 = 0.7821. Previously was
%0.6
%f_sw_calc = 0.8879/0.9834; %test to match the 2nd period better for Hist-Aer

tr_atmos_constant_val = 0.95;
tr_atmos_constant_val = 0.7821;
tr_atmos_constant_val = 0.634; %0.76 used in Seinfeld and Pandis (Section 24.8.2), so not too far off that
tr_atmos_constant_val = 0.76;
tr_atmos_constant_val = 0.8; %0.83
tr_atmos_constant_val = 0.83; %0.83
tr_atmos_constant_val = 0.89; 

%tr_atmos_constant_val = 0.734; %attempt to get this to match using cloudy-sky trans scaling only.

cf_min = 0.01;
i_Liu = 1;


fac_tr = 1/1.1;
fac_tr = 1;

fac_SWin=1.7;
fac_SWin=1;

fac_Nd = 1;
%fac_Nd = 2;%1.3;

fac_cf = 1;
%fac_cf = 1.2;

fac_SW_calc=1;
%fac_SW_calc=1.4;

%% Atmospheric transmissivity
switch tr_for_SW_bias
    case 'constant value'
        tr_atmos_mean = tr_atmos_constant_val; %if constant value
    case 'time mean from model'
        SW_thresh=200;
        %transmission_atmos = SW_clean_clear_PD_ALL ./ SW_in_PD_ALL;
        %transmission_atmos(SW_in_PD_ALL<SW_thresh)=0.5;
        %tr_atmos_mean = meanNoNan(transmission_atmos,3);
        
    case {'time varying monthly','time varying annual mean'}
        transmission_atmos_up = sw_TOA_up_clear_dat_ens_mean ./ sw_surf_up_clear_dat_ens_mean;
        %Upwelling from at TOA divided by upwelling at surface in
        %clear-skies.
        
        transmission_atmos_down = sw_surf_down_clear_dat_ens_mean ./ SW_down_TOA_dat_ens_mean;
        %Downwelling to the surface divided by downwelling at TOA in
        %clear-skies.
        
        %So, will get two measures of transmissivity here. The TOA
        %calculation uses a cloud top SW downwelling, so could just use
        %the surface downwelling for this (although may underestimate a
        %little).
        %Then can provide the upwelling trans value to the SW calc
        
        %Actually, can only really use the downwards transmission since upwards includes
        %the scattered light from the sun and not just that originating
        %from the surface.
end

%% Surface albedo
switch albedo_type
    case 'constant'
        % Also need the clear-sky albedo - can estimate from cloud free regions?
        surf_albedo = 0.15; %Guess for now - value for CF quite sensitive to this (ranging up to 30%), so should prob
        %check what this is from the model
        surf_albedo = 0.12; %value estimated from the 12 noon snapshot from VOCALS reginal run - N.B. - this uses the estimated transmission since
        %it is the surface albedo assuming no atmosphere above
        %surf_albedo = 0.0;
        surf_albedo = 0.09;
        
    case 'calculated monthly'
        %Estimate of surface albedo using surface fluxes in clear-sky
        surf_albedo  = dat_PI.SW_up_surf_clear_sky ./ dat_PI.SW_down_surf_clear_sky;
        
    case 'calculated annual mean'
        surf_albedo_monthly  = sw_surf_up_clear_dat_ens_mean ./ sw_surf_down_clear_dat_ens_mean;
        
        clear surf_albedo_ann
        for iy=1:length(years_sw_calc)
            istart=(iy-1)*12+1;
            %surf_albedo_ann(iy,:,:) = meanNoNan(surf_albedo_monthly(istart:istart+11,:,:),1,'mean',1,1,SW_down_TOA_dat_ens_mean(istart:istart+11,:,:));
            surf_albedo_ann(iy,:,:) = min(surf_albedo_monthly(istart:istart+11,:,:),[],1);
        end
        %trans_orig = repmat(trans_orig_ann,[1 1 1 12]);
        surf_albedo = NaN*ones(size(surf_albedo_monthly));
        for im=1:size(surf_albedo_monthly,1)
            iy = ceil(im/12);
            surf_albedo(im,:,:) = surf_albedo_ann(iy,:,:);
        end
        
        
        
end


