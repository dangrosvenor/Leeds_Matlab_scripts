% SW_calcs_Hawaii_setup
% Called from SW_calcs_Hawaii_RUN.m

%% 2) Set up stuff for the SW TOA calcs from monthly means 

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
albedo_type = 'calculated all time steps';
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
