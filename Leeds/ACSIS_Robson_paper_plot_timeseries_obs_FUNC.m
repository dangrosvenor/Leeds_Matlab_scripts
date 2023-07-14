
%For some variables there is a choice of obs - the line below sets the one
%required.
eval(['obs_str_in = obs_str' obs_var_end_str ';']);


switch var_ukesm   
        
    case 'SO2_low_anthropogenic_emissions'
                % DEEP-C SW dataset for SW
        obs_str_out='DEEP-C';
        years_obs_out = [1985:2014];
        
        dat_obs_out = dat_deepc.sw_up_toa;
        
        lat_obs_out = gcm_Plat2D_UM;
        lon_obs_out = gcm_Plon2D_UM;
        gcm_area_obs_out = gcm_area_UM;
        
        %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';  
        yr_start_trend_obs_out = years_obs_out(1); 
        yr_end_trend_obs_out = 2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend); 
        
        clear opts; opts.screen_type_optional = land_ocean;
        [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
        
%         if iscreen_land==1
%             [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs_out,lat_obs_out,lon_obs_out);
%         else
%             dat_tmp = dat_obs_out;
%         end
        
        fscale=1e13; fscale_str=['\times10^{-13}'];
        
    case {'SW_up_TOA'}
        % DEEP-C SW dataset for SW
        obs_str_out='DEEP-C';
        years_obs_out = [1985:2014];
        
        dat_obs_out = dat_deepc.sw_up_toa;
        
        lat_obs_out = gcm_Plat2D_UM;
        lon_obs_out = gcm_Plon2D_UM;
        gcm_area_obs_out = gcm_area_UM;
        
        %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';  
        yr_start_trend_obs_out = years_obs_out(1); 
        yr_end_trend_obs_out = 2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);  
        
        clear opts; opts.screen_type_optional = land_ocean;
        [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
        
        %if iscreen_land==1
        %    [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs_out,lat_obs_out,lon_obs_out);
        %else
        %    dat_tmp = dat_obs_out;        
        %end
        
        fscale=1e2; fscale_str=['\times10^{-2}'];

        
    case {'rsut'}
        % DEEP-C SW dataset for SW
        obs_str_out='DEEP-C';
        years_obs_out = [1985:2014];
        
        dat_obs_out = dat_deepc.sw_up_toa;
        
        lat_obs_out = gcm_Plat2D_UM;
        lon_obs_out = gcm_Plon2D_UM;
        gcm_area_obs_out = gcm_area_UM;
        
        %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';  
        yr_start_trend_obs_out = years_obs_out(1); 
        yr_end_trend_obs_out = 2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);  
        
        clear opts; opts.screen_type_optional = land_ocean;
        [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
        
%         if iscreen_land==1            
%             [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs_out,lat_obs_out,lon_obs_out);
%         else
%             dat_tmp = dat_obs_out;  
%         end
        
        fscale=1e2; fscale_str=['\times10^{-2}'];  
        
    case {'SWTOA Calc'}
        % DEEP-C SW dataset for SW
        obs_str_out='Calculated';
        years_obs_out = dat_ukesm_swcalc.years_ukesm_1d;
        
        dat_obs_out = dat_ukesm_swcalc.dat_monthly; %Need monthly data here.
        
        lat_obs_out = gcm_Plat2D_UM;
        lon_obs_out = gcm_Plon2D_UM;
        gcm_area_obs_out = gcm_area_UM;
        
        %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';  
        yr_start_trend_obs_out = yr_start_trend_box(iconstant_trend); %years_obs_out(1); 
        yr_end_trend_obs_out = yr_end_trend_box(iconstant_trend); %2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend); 
        
        yr_start_trend_obs_out = yr_start_trend_box; %years_obs_out(1); 
        yr_end_trend_obs_out = yr_end_trend_box; %
        
        
        clear opts; opts.screen_type_optional = land_ocean;
        [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
        
%         if iscreen_land==1            
%             [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs_out,lat_obs_out,lon_obs_out);
%         else
%             dat_tmp = dat_obs_out;  
%         end
        
        fscale=1e2; fscale_str=['\times10^{-2}'];        
        
    case {'calipso_total_cloud_amount','clt'}
        
        switch obs_str_in %
            case 'CCI'
                % CCI total CF dataset
                obs_str_out='ESA CCI';
                
                years_obs_out = years_obs_CCI;
                dat_obs_out = cci_dat.cfc; %need monthly values here
                
                lat_obs_out = gcm_Plat2D_CCI;
                lon_obs_out = gcm_Plon2D_CCI;
                gcm_area_obs_out = gcm_area_CCI;
                
                yr_start_trend_obs_out = 1985;
                yr_end_trend_obs_out = 2014;
                
                lmask_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ESA_CCI_landmask_data_cfc.mat';
                load(lmask_file,'lmask_out');
                
                clear opts; opts.screen_type_optional = land_ocean; opts.lmask_in_optional = lmask_out;
                [dat_tmp,lmask_out2] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
        
            case 'PATMOSx'
                % PATMOS Norris anomaly dataset
                obs_str_out='PATMOSx';
                patmos_baseline=0.66; %Norris data is monthly anomalies relative to a global spatial and time-mean CF
                 %(I think).
                
                years_obs_out = years_obs_PATMOS;
                dat_obs_out = patmos_baseline + patmos_dat.totcf_anom_PATMOS/100; %need monthly values here
                
                lat_obs_out = gcm_Plat2D_UM;
                lon_obs_out = gcm_Plon2D_UM;
                gcm_area_obs_out = gcm_area_UM;
                
                yr_start_trend_obs_out = 1983;
                yr_end_trend_obs_out = 2009;
                
                %Grid is the same as for UM
                clear opts; opts.screen_type_optional = land_ocean;
                [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
                
            case 'ISCCP'
                % ISCCP Norris anomaly dataset
                obs_str_out='ISCCP';
                isccp_baseline=0.66; %Norris data is monthly anomalies relative to a global spatial and time-mean CF
                 %(I think).
                
                years_obs_out = years_obs_ISCCP;
                dat_obs_out = isccp_baseline + isccp_dat.totcf_anom_ISCCP/100; %need monthly values here
                
                lat_obs_out = gcm_Plat2D_UM;
                lon_obs_out = gcm_Plon2D_UM;
                gcm_area_obs_out = gcm_area_UM;
                
                yr_start_trend_obs_out = 1983;
                yr_end_trend_obs_out = 2009;
                
                %Grid is the same as for UM
                clear opts; opts.screen_type_optional = land_ocean;
                [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);                              
                
        end
                
        
        fscale=1e4; fscale_str=['\times10^{-4}'];
        
    case 'lwp'         
                obs_str_out='MAC';
                                
                years_obs_out = years_obs_MAC;
                dat_obs_out = MAC_monthly; %need monthly values here
                
                lat_obs_out = gcm_Plat2D_UM;
                lon_obs_out = gcm_Plon2D_UM;
                gcm_area_obs_out = gcm_area_UM;
                
                yr_start_trend_obs_out = 1988;
                yr_end_trend_obs_out = 2014;
                
                %Grid is the same as for UM
                clear opts; opts.screen_type_optional = land_ocean;
                [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
        
case {'Nd_cf_weighted_UKESM','scldncl','Nd_clw_weighted_ESGF','Nd_clw_weighted_ESGF_no_dz'}
    iplot_obs01=0;
    iplot_obs02=0; 
    
    %MODIS dataset
    switch Nd_obs
        case {'MODIS','both'}
            obs_str_out = 'MODIS';
            years_obs_out = years_MODIS2;
            
            dat_obs_out = Nd_monthly_MODIS;
            
            lat_obs_out = gcm_Plat2D_AMSRE;
            lon_obs_out = gcm_Plon2D_AMSRE;
            gcm_area_obs_out = gcm_area_AMSRE;
            
            yr_start_trend_obs_out = years_MODIS2(1);
            yr_end_trend_obs_out = years_MODIS2(end);            
            
            clear opts; opts.screen_type_optional = land_ocean;
            [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);                        
            
                 
            %fscale=1e4; fscale_str=['\times10^{-4}'];
            fscale=1; fscale_str=[''];       
            
            iplot_obs02=1; %flag to plot MDOIS
            
            
            
    end
    
    
    switch Nd_obs
        case {'CCI','both'}
            
            % CCI Nd dataset - also add the MODIS option here.
            obs_str='ESA CCI';
            %years_obs = unique(cci_dat.Y_out);
            years_obs = unique(dat_Nd_CCI.Y_out);
            
            
            
            dat_obs = dat_Nd_CCI.N_T268_cf20_cot5_re3; %
            
            lat_obs = gcm_Plat2D_CCI;
            lon_obs = gcm_Plon2D_CCI;
            gcm_area_obs = gcm_area_CCI;
            
            yr_start_trend_obs = 1985;
            yr_end_trend_obs = 2014;
            
            lmask_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ESA_CCI_landmask_data_Nd.mat';
            load(lmask_file,'lmask_out');
            lmask_cci = repmat(lmask_out,[1 1 size(dat_obs,1)]);
            lmask_cci = permute(lmask_cci,[3 1 2]);
            
            
            %fac = 90/50;
            fac=1;
            %offset=40; %seems to work best to bring it inline with MODIS
            offset=0;
            dat_tmp_cci = fac * dat_obs .* lmask_cci + offset;
            
            [obs_monthly_box,obs_annual_box,obs_annual_map] ... 
    = ACSIS_Robson_paper_process_obs(dat_tmp_cci,lat_obs,lon_obs,LAT_val,LON_val,years_obs,gcm_area_obs,season);
       
            
            %fscale=1e4; fscale_str=['\times10^{-4}'];
            fscale=1; fscale_str=[''];
            
            inc_obs_01=1;            
            iplot_obs01=1; %flag to plot CCI
           
            
    end
    
    case {'ts'}
        switch obs_str_in
            case 'HadCRUTv5'
                obs_str_out='HadCRUTv5';
                
                years_obs_out = years_obs_ts;
                dat_obs_out = 273.15 + permute(ts_dat.tot_tas_anom_HadCRUT,[3 1 2]); %need monthly values here
                
                lat_obs_out = gcm_Plat2D_UM;
                lon_obs_out = gcm_Plon2D_UM;
                gcm_area_obs_out = gcm_area_UM;
                
                yr_start_trend_obs_out = 1850;
                yr_end_trend_obs_out = 2021;
                
                %Grid is the same as for UM
                clear opts; opts.screen_type_optional = land_ocean;
                [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
                
            case {'AMIP','HadISST'} %AMIP is run with HadISST temperatures.
                % surface temperatures from AMIP run - change to actual HADISST
                % data
                %obs_str_out = 'UKESM1 atmos-only';
                obs_str_out = 'UKESM1-AMIP';
                %years of the data itself
                years_obs_out = [1979:2014];
                
                dat_obs_out = ts_dat.dat_ens_mean;
                
                lat_obs_out = ts_dat.gcm_Plat2D_UM;
                lon_obs_out = ts_dat.gcm_Plon2D_UM;
                gcm_area_obs_out = gcm_area_UM; %same as model if using AMIP.
                
                %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';
                yr_start_trend_obs_out = 1985; %ts_dat.years_ukesm_1d(1);
                yr_end_trend_obs_out = 2014; %ts_dat.years_ukesm_1d(end); %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
                
                clear opts; opts.screen_type_optional = land_ocean;
                [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
                
                %         if iscreen_land==1
                %             [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs_out,lat_obs_out,lon_obs_out);
                %         else
                %             dat_tmp = dat_obs_out;
                %         end
                
                %fscale=1; fscale_str=['\times10^{-2}'];
                
        end
        
        
  case {'od550aer','od550tot'}
        % MODIS AOD Aqua C6.1 2001-2021
        obs_str_out='MODIS Aqua';
        years_obs_out = years_obs_MODIS; %[2001:2021];
        
        dat_obs_out = modis_aod.aod;
        
        lat_obs_out = gcm_Plat2D_MODIS;
        lon_obs_out = gcm_Plon2D_MODIS;
        gcm_area_obs_out = gcm_area_MODIS;
        
        %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';  
        yr_start_trend_obs_out = years_obs_out(1); 
        yr_end_trend_obs_out = 2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);  
        
        clear opts; opts.screen_type_optional = land_ocean;
        [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_obs_out,lat_obs_out,lon_obs_out,opts);
        
        %if iscreen_land==1
        %    [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs_out,lat_obs_out,lon_obs_out);
        %else
        %    dat_tmp = dat_obs_out;        
        %end
        
        fscale=1; fscale_str=[''];       
            
        
end

[obs_monthly_box_out,obs_annual_box_out,obs_annual_map_out] ... 
    = ACSIS_Robson_paper_process_obs(dat_tmp,lat_obs_out,lon_obs_out,LAT_val,LON_val,years_obs_out,gcm_area_obs_out,season);


var_list={'obs_str' 
    'years_obs'
    'dat_obs'
    'lat_obs'
    'lon_obs'
    'gcm_area_obs'
    'yr_start_trend_obs'
    'yr_end_trend_obs'
    'obs_monthly_box'
    'obs_annual_box'
    'obs_annual_map'};


for ivar=1:length(var_list)
   eval_str = [var_list{ivar} obs_var_end_str ' = ' var_list{ivar} '_out;'];
   eval(eval_str);
end