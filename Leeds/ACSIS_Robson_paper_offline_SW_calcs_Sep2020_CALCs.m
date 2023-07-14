%This script calculates the SWTOA estimates for all vars varying and for
%each one kept constant (for every month globally). Then saves in a .mat file - one for each period
%(since uses different constant values for the different periods).
%Runs over a loop to do the two periods
%Then run
%   ACSIS_Robson_paper_offline_SW_calcs_Sep2020.m
%to plot and save
%the trends for the regional box (annual means) and make the table.
% Don't need to run ACSIS_Robson_paper_load_data.m first

%----------------------------------------------------------------------
%----- SW calculation/recreation using model monthly means
%----------------------------------------------------------------------

%Need to select which trend period we want to calculate and save for
%Need to run for each period requried.

%isingle_ens=0;
%isingle_ens=1; ens_inds=2; %if want to load individual ensemble members - choose which one here.

%MIP = 'CMIP'; expt_multi = {''}; ens_inds=[1:16]; isingle_ens=1; %(UKESM)
MIP = 'DAMIP'; expt_multi = {'DAMIP_hist-aer','DAMIP_hist-GHG','DAMIP_hist-nat'}; isingle_ens=1; ens_inds=[1:4]; %now loops through and does a running average of all the values as well as saving individual .mat files for each member
MIP = 'AerChemMIP'; expt_multi = {'UKESM1-AerChemMIP_control'}; isingle_ens=1; ens_inds=[1:3]; %now loops through and does a running average of all the values as well as saving individual .mat files for each member
MIP = 'AerChemMIP'; expt_multi = {'AerChemMIP_hist-piAer'}; isingle_ens=1; ens_inds=[1:3]; %now loops through and does a running average of all the values as well as saving individual .mat files for each member
%MIP = 'AerChemMIP'; expt_multi = {'AerChemMIP_hist-AerProxy'}; isingle_ens=1; ens_inds=[1:3]; %now loops through and does a running average of all the values as well as saving individual .mat files for each member
%MIP = 'HADGEM3_GC31_LL'; expt_multi = {'HADGEM3_GC31_LL'}; ens_inds=[1:4]; isingle_ens=1; %(HADGEM)

Nens = length(ens_inds);

yr_start_trend_box = [1870 1985]; yr_end_trend_box = [1970 2014]; %
yr_start_trend_box = [1850 1971]; yr_end_trend_box = [1970 2014]; %
%which trends to plot is chosen later - search for iplot_trend

for iperiod=1:2
    
    if iperiod==1
        iconstant_trend=1; iplot_trend=[1 0]; yr_start_plot = 1850;
    else
        iconstant_trend=2; iplot_trend=[0 1]; yr_start_plot = 1979;
        %iconstant_trend=1; iplot_trend=[1 1]; yr_start_plot = 1850; yr_end_plot = 2018;
    end
    
    itrend_lab = iconstant_trend;
    
    %iadd_trend_str=1; %whether to add the trend value to the legend
    iadd_trend_str=0;
    ipad_legend = 0; %Pad the first entry to avoid box cutting into the text when have
    %superscripts
    
    isave_plot_SWTOA_calcs_draft = 0;
    years_sw_calc = [1984:2014];
    years_sw_calc = [1850:2014];
    
    
    
    
    %% Loop over various models
    
    
    %Loop over all models
    
    for idamip_run=1:length(expt_multi)
        
        iens_count = 0;
        for iens=ens_inds
            iens_count = iens_count + 1;
            
            expt_str = expt_multi{idamip_run};
            
            %% 1) Load ensemble monthly mean data
            
            
            
            %SW down TOA
            %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_SW_down_TOA_partial_SW_down_TOA.mat';
            load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsdt.mat';
            
            load(load_file,'gcm_Plat2D_UM','gcm_Plon2D_UM','gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM');
            
            %Don't have this for the DAMIP runs yet - but can probably assume that
            %this is the same as in the UKESM since should not be affected by the
            %planet. Worth checking, though.
            
            %Pick out just the years for period 1 from others (all months for 1850-2014)
            %(1984-1850)*12 = 1608
            load(load_file,'years_ukesm_1d');
            iy01 = find(years_ukesm_1d==years_sw_calc(1));
            iy02 = find(years_ukesm_1d==years_sw_calc(end));
            %replicate the years for each month
            clear years_sw_calc_monthly months_sw_calc_monthly
            for iy=1:length(years_sw_calc)
                istart = (iy-1)*12+1;
                iend=istart+11;
                years_sw_calc_monthly(istart:iend)=years_sw_calc(iy);
                months_sw_calc_monthly(istart:iend)=[1:12];
            end
            %t_inds = 1609:1980;
            t_inds = (iy01-1)*12+1:(iy02-1)*12+12; %data is monthly, so convert indices.
            t_inds_annual = iy01:iy02; %data is monthly, so convert indices.
            
            %find the indices for the trend period to keep constant
            iy01_trend = find(years_sw_calc_monthly==yr_start_trend_box(iconstant_trend));
            iy02_trend = find(years_sw_calc_monthly==yr_end_trend_box(iconstant_trend));
            inds_trend = iy01_trend:iy02_trend;
            
            %Load ensemble mean data.
            
            %To avoid loading the whole dat_ens array can do :-
            % mat_obj = matfile(load_file)
            % and then :- dat = mat_obj.dat_ens(iens,t_inds,:,:)
            % Also, to determine size of an array :- size(mat_obj,'dat_ens')
            
            mat_obj = matfile(load_file);
            if isingle_ens==1
                SW_down_TOA_dat_ens_mean = squeeze(mat_obj.dat_ens(iens,t_inds,:,:));
            else
                SW_down_TOA_dat_ens_mean = mat_obj.dat_ens_mean(t_inds,:,:);
            end
            
            % --- SW TOA upwelling CLEAR-SKY - for calculation of transmissivity ---
            var_DAMIP = 'rsutcs';
            switch MIP
                case {'DAMIP','AerChemMIP','HADGEM3_GC31_LL'}
                    fscale = 1;
                    %script to set the load file - checks for special case
                    %of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP                                        
                    
                otherwise
                    fscale = 1;
                    model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
            end
            
            %run the function to calculate the aerosol proxy for AerChemMIP
            [sw_TOA_up_clear_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);            
                        
            
            clear sw_TOA_up_clear_dat_annual
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                sw_TOA_up_clear_dat_annual(iy,:,:) = meanNoNan(sw_TOA_up_clear_dat_ens_mean(istart:istart+11,:,:),1);                
            end
            
            
            % --- SW surface downwelling CLEAR-SKY - for calculation of transmissivity and surface albedo ---
            var_DAMIP = 'rsdscs';
            switch MIP
                case {'DAMIP','AerChemMIP','HADGEM3_GC31_LL'}
                    fscale = 1;
                    %script to set the load file - checks for special case
                    %of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP 
                    
                    
                otherwise
                    fscale = 1;
                    model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
            end
            
            %run the function to calculate the aerosol proxy for AerChemMIP
            [sw_surf_down_clear_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);                                                    
            
            
            % --- SW surface upwelling CLEAR-SKY - for calculation of transmissivity and surface albedo ---
            var_DAMIP = 'rsuscs';
            switch MIP
                case {'DAMIP','AerChemMIP','HADGEM3_GC31_LL'}
                    fscale = 1;
                    %script to set the load file - checks for special case
                    %of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP 
                    
                otherwise
                    fscale = 1;
                    model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
            end
            
            %run the function to calculate the aerosol proxy for AerChemMIP
            [sw_surf_up_clear_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);                                                    
                                   
            
            
            % Nd
            switch MIP
                case {'DAMIP'}
                    var_DAMIP = 'Nd_clw_weighted_ESGF';
                    fscale = 1e-6; %convert from m^-3 to cm^-3
                    
                    %script to set the load file - checks for special case
                    %of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP 
                    
                    
                    
                case {'AerChemMIP'}
                    var_DAMIP = 'Nd_clw_weighted_ESGF_total_column_to_zdomain_top';
                    fscale = 1e-6; %convert from m^-3 to cm^-3
                    
                   %script to set the load file - checks for special case of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP 
                      
                otherwise
                    %var_DAMIP = var_ukesm;
                    %fscale_DAMIP = dat_ukesm.fscale;
                    
                    fscale = 1;
                    
                    %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';
                    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_scldncl.mat';
            end  
            
            %run the function to calculate the aerosol proxy for AerChemMIP
            [Nd_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);                                                    
            
            
           
            % --- total CF ---
            switch MIP
                case {'DAMIP','AerChemMIP','HADGEM3_GC31_LL'}
                    var_DAMIP = 'clt';
                    fscale = 0.01;                    
                    %script to set the load file - checks for special case of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP                     
                otherwise
                    fscale = 0.01;
                    
                    %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat';
                    %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_total_cloud_amount.mat';
                    model_str = 'UKESM1'; load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_clt.mat';
            end
            
            %run the function to calculate the aerosol proxy for AerChemMIP
            [totcf_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);                                                    
                                   
            
            %low CF - don't think model low CF this exists for ESGF - do have CALIPSO
            %low CF, though - in /badc/cmip6/data/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical
            %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat';
            %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat';
            %load(load_file,'dat_ens_mean'); %size [1980         144         192];
            %1980 is all 12 months of all 165 years.
            %lowcf_dat_ens_mean = dat_ens_mean(t_inds,:,:);
            
            
            
            % --- LWP (kg/m2) ---
            %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_LWP_2-391.mat';
            
            %1980 is all 12 months of all 165 years.
            
            fscale = 1e-3; %convert to g/m2 from kg/m2
            var_DAMIP = 'lwpic';
            switch MIP
                case {'DAMIP','AerChemMIP','HADGEM3_GC31_LL'}
                    %fscale = 0.01;                    
                    %script to set the load file - checks for special case of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP                                         
                otherwise
                    %fscale = 0.01;
                    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_lwpic.mat';
                    
            end
            
            %run the function to calculate the aerosol proxy for AerChemMIP
            [lwp_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);                                                    
                                   
            
        % --- Aerosol optical depth
            var_DAMIP = 'od550aer';            
            switch MIP
                case {'DAMIP'}
                    fscale = 1;                    
                    %script to set the load file - checks for special case of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP                                         
                case {'AerChemMIP','HADGEM3_GC31_LL'}    
                    load_file = 'dummy'; %Not actually needed if using clear-sky SW.
                case 'CMIP'                
                    fscale = 1;
                    model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
                otherwise
                    load_file = 'dummy';
            end                       
            
            switch load_file
                case 'dummy'
                    aod_dat_ens_mean = ones(size(lwp_dat_ens_mean))*NaN; %Not actually needed if using clear-sky SW.
                otherwise
                    %run the function to calculate the aerosol proxy for AerChemMIP
                    [aod_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);
                    
            end
            
       % --- Absorbing aerosol optical depth
            var_DAMIP = 'abs550aer';
            
            switch MIP
                case {'DAMIP'}
                    fscale = 1;
                    %script to set the load file - checks for special case of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP                                                             
                case {'AerChemMIP','HADGEM3_GC31_LL'}
                    load_file = 'dummy';
                case 'CMIP'
                    fscale = 1;
                    model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
                otherwise
                    load_file = 'dummy'; %default
            end
            
            switch load_file
                case 'dummy'
                    aaod_dat_ens_mean = ones(size(lwp_dat_ens_mean))*NaN;
                otherwise
                    %run the function to calculate the aerosol proxy for AerChemMIP
                    [aaod_dat_ens_mean] = ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);                                        
            end
            
         % --- Dust optical depth
%             var_DAMIP = 'od550dust';
%             switch MIP
%                 case {'DAMIP','AerChemMIP','HADGEM3_GC31_LL'}
%                     fscale = 1;
%                     load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str '_all_' var_DAMIP '.mat'];
%                     
%                 otherwise
%                     fscale = 1;
%                     model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
%             end
%             
%             mat_obj = matfile(load_file);
%             if isingle_ens==1
%                 dust_aod_dat_ens_mean = squeeze(fscale *  mat_obj.dat_ens(iens,t_inds,:,:));
%             else
%                 dust_aod_dat_ens_mean = fscale * mat_obj.dat_ens_mean(t_inds,:,:); %convert from % to 0 to 1 form.
%             end            
           
            
            % -- SW TOA up (to compare against the calculation) --
            fscale = 1;
            var_DAMIP = 'rsut';
            switch MIP
                case {'DAMIP','AerChemMIP','HADGEM3_GC31_LL'}
                    %fscale = 0.01;                    
                    %script to set the load file - checks for special case of aerosol only proxy using control minus piAer
                    ACSIS_Robson_paper_offline_SW_set_loadfile_AerChemMIP                                         
                otherwise
                    %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat';
                    %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_SW_up_TOA.mat';
                    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_rsut.mat';
            end
            
            %run the function to calculate the aerosol proxy for AerChemMIP
            [SW_up_TOA_dat_ens_monthly,SW_up_TOA_dat_annual_ens_mean,SW_up_TOA_dat_ens_mean,SW_up_TOA_dat_annual_ens] =...
                ACSIS_Robson_paper_offline_SW_calc_AerProxy_AerChemMIP(load_file,load_file2,isingle_ens,fscale,iens,t_inds,t_inds_annual,expt_str,var_DAMIP);
                                                                      
            %
%             mat_obj = matfile(load_file);
%             if isingle_ens==1
%                 SW_up_TOA_dat_annual_ens_mean = squeeze( mat_obj.dat_annual_ens(iens,t_inds_annual,:,:) );
%                 SW_up_TOA_dat_ens_monthly = squeeze(mat_obj.dat_ens(iens,t_inds,:,:));
%             else
%                 %load(load_file,'dat_ens_mean','dat_annual_ens','dat_annual'); %size [1980         144         192];
%                 %1980 is all 12 months of all 165 years.
%                 SW_up_TOA_dat_annual_ens_mean = mat_obj.dat_annual(t_inds_annual,:,:);
%                 SW_up_TOA_dat_ens_monthly = squeeze(mat_obj.dat_ens(t_inds,:,:));
%             end
%             
%             SW_up_TOA_dat_ens_mean = mat_obj.dat_ens_mean(t_inds,:,:);
%             SW_up_TOA_dat_annual_ens = mat_obj.dat_annual_ens(:,t_inds_annual,:,:);
            
            
            %% 2) Set up stuff for the SW TOA calcs from monthly means - likely inaccuracies due to using monthly mean SWTOA down, etc.
            %N.B. - can get the clear-sky downwelling at the surface to calc this if
            %needed.
            
            %f_upscatter = 0.4;
            
            trans_clear_sky = 0.95;
            trans_clear_sky = 1.0;
            
            f_up_method = 'sza detailed';
            
            switch f_up_method
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
                    
                    %Calculate SZA just using the 15th of each month for one year
                    %at 0 deg lon.
                    
                    dt=0.5;
                    savedir = '/home/disk/eos15/d.grosvenor/UM/ACSIS_Nd_trends/';
                    savename_sza = [savedir 'sza_' num2str(dt) '_SW_calcs.mat'];
                    
                    i_calc = 0;
                    if i_calc==1
                        
                        %mat_times_monthly = datenum(years_sw_calc_monthly,months_sw_calc_monthly,15,time_LT_hr,time_LT_min,0);
                        
                        HH = [0:dt:24-dt]; %calculate every 30 mins
                        HH2 = repmat(HH,[12 1]); %Replicate to every month
                        DD2 = 15*ones(size(HH2)); %15th day of  month
                        MM2 = repmat([1:12],[length(HH) 1]); MM2 = MM2';
                        YY2 = 1850*ones(size(HH2)); %year 1850
                        %mat_times_monthly = datenum(1850*ones([1 12]),[1:12],15,time_LT_hr,time_LT_min,0);
                        %mat_times_ALL = repmat(mat_times_monthly',[1 size(SW_up_TOA_dat_ens_mean,2)]);
                        
                        mat_times = datenum(YY2,MM2,DD2,HH2,0,0);
                        lats = gcm_Plat2D_UM(:,1);
                        mat_times_ALL = repmat(mat_times,[1 1 length(lats)]);
                        lats_ALL = repmat(lats,[1 size(mat_times)]); lats_ALL = permute(lats_ALL,[2 3 1]);
                        lons_ALL = zeros(size(lats_ALL(:))); %Use UTC, so specify 0 deg lon
                        sza = sun_pos(mat_times_ALL,lats_ALL,lons_ALL); sza(sza>90) = 90;
                        %sza_ALL = repmat(sza,[size(SW_up_TOA_dat_ens_mean,1)/12 1 size(SW_up_TOA_dat_ens_mean,3)]);
                        
                        %[12 Nhrs 144]
                        me_cos = meanNoNan(cos(pi/180*sza),2);
                        %me_sza = 180/pi * acos(me_cos);
                        
                        %replicate to all years and lons
                        me_cos_ALL = repmat(me_cos,[size(SW_up_TOA_dat_ens_mean,1)/12 1 size(SW_up_TOA_dat_ens_mean,3)]);
                        
                        
                        save(savename_sza,'sza','me_cos','me_cos_ALL','dt','-V7.3');
                        
                    else
                        load(savename_sza);
                    end
                    
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
            
            %Value used in draft for DAMIP, etc. - shoudl prob keep same
            %value for all if possible?
            tr_atmos_constant_val = 0.734; %attempt to get this to match using cloudy-sky trans scaling only.            
            tr_atmos_constant_val = 0.726;
            
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
                    surf_albedo  = sw_surf_up_clear_dat_ens_mean ./ sw_surf_down_clear_dat_ens_mean;
                    
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
            
            
            
            %% Do actual calc
            time_mean=0;
            time_mean=1;
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            
            lwp_type = 'new';
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %already in in-cloud LWP
            end
            
            nd = Nd_dat_ens_mean;
            
            switch tr_for_SW_bias
                case 'constant value'
                    trans_orig = tr_atmos_constant_val; %if constant value
                    %sw = meanNoNan(SW_in,3);
                   % sw = fac_SWin * SW_down_TOA_dat_ens_mean .* trans_orig;
                    sw = fac_SWin * SW_down_TOA_dat_ens_mean; 
                    
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
                        
            
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = sw_TOA_up_clear_dat_ens_mean;
            %[Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model,SWTOA_clear_sky_calc_model] = calc_SW2(cf,lwp*1e3,nd,sw,surf_albedo,trans_orig,cf_min,0,NaN,i_Liu);
            [Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model,SWTOA_clear_sky_calc_model] = calc_SW3(opts,aod_dat_ens_mean,aaod_dat_ens_mean,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,surf_albedo,trans_orig);            

            SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual SWTOA_clear_sky_calc_model_annual
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual(iy,:,:) = meanNoNan(SWTOA_calc_model(istart:istart+11,:,:),1);
                SWTOA_clear_sky_calc_model_annual(iy,:,:) = meanNoNan(SWTOA_clear_sky_calc_model(istart:istart+11,:,:),1);
            end
            
            clear surf_albedo_annual
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                surf_albedo_annual(iy,:,:) = meanNoNan(surf_albedo(istart:istart+11,:,:),1);
            end
            
            
            %% Do some maps comparing calculated to actual
            i_plot_maps=0;
            if i_plot_maps==1
                
                
                
                irestrict_domain_DRIVER=1;
                ioverride_proj_type=1;
                proj_type_DRIVER='ortho';
                ioverride_LAT_plots=0;
                iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
                ioverride_ticks_DRIVER=1;
                icoarse_grain=0;
                icontour_DRIVER=0; cont_col_str_DRIVER='';
                time_round=0; time_format_str='';
                isave_plot=0;
                iplot_wind_arrows=0;
                
                iyear = find(years_sw_calc==1965);
                iy=iyear-5:iyear+5;
                %iy=156:165;
                iy=1:11;
                
                    ioverride_LAT_plots=0;
                    dat_modis = meanNoNan(SWTOA_calc_model_annual(iy,:,:),1); 
                    %var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' iy=' num2str(iy)];
                    var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
                    %run plotting script
                    figure
                    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
                    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
                    caxis([0 250]);
                
                    ioverride_LAT_plots=0;
                    dat_modis = meanNoNan(SW_up_TOA_dat_annual_ens_mean(iy,:,:),1); 
                    var_UM = 'Actual SWTOA'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
                    %run plotting script
                    figure
                    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
                    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
                    caxis([0 250]);

                    
                    ioverride_LAT_plots=0;
                    dat_modis = meanNoNan(SWTOA_clear_sky_calc_model_annual(iy,:,:),1); 
                    var_UM = 'Calculated SWTOA clear-sky'; subtitle_str = [var_UM 'trans_clear_sky=' num2str(trans_clear_sky) ' m=' num2str(m) ' c=' num2str(c)];
                    %run plotting script
                    figure
                    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
                    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
                    caxis([0 75]);
                
                    ioverride_LAT_plots=0;
                    dat_modis = meanNoNan(sw_TOA_up_clear_dat_annual(iy,:,:),1); 
                    var_UM = 'Actual SWTOA clear-sky'; subtitle_str = [var_UM 'trans_clear_sky=' num2str(trans_clear_sky) ' m=' num2str(m) ' c=' num2str(c)];
                    %run plotting script
                    figure
                    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
                    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
                    caxis([0 75]);  
                    
                    ioverride_LAT_plots=0;
                    dat_modis = meanNoNan(SWTOA_clear_sky_calc_model_annual(iy,:,:),1) - meanNoNan(sw_TOA_up_clear_dat_annual(iy,:,:),1); 
                    var_UM = 'Bias SWTOA clear-sky'; subtitle_str = [var_UM 'trans_clear_sky=' num2str(trans_clear_sky) ' m=' num2str(m) ' c=' num2str(c)];                    
                    %run plotting script
                    figure
                    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
                    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
                    caxis([- 20 20]); 
                
%                 ioverride_LAT_plots=0;
%                 dat_modis = meanNoNan(surf_albedo_annual(iy,:,:),1); var_UM = 'Annual mean surface albedo'; subtitle_str = [var_UM ' iy=' num2str(iy)];
%                 %run plotting script
%                 figure
%                 UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                 lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                 caxis([0 0.2]);
                
%                 ioverride_LAT_plots=0;
%                 dat_modis = meanNoNan(trans_orig_ann(iy,:,:),1); var_UM = 'Annual mean transmissivity'; subtitle_str = [var_UM ' iy=' num2str(iy)];
%                 %run plotting script
%                 figure
%                 UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                 lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                 caxis([0.6 0.9]);
                
                %montly plots
                %
                
                ifind = find(years_sw_calc_monthly==1970);
                ims=ifind:ifind+11;
                
%                 
%                 for im=ims
%                     
%                     imon = im-ims(1)+1
%                     %
%                     %         ioverride_LAT_plots=0;
%                     %         dat_modis = squeeze(SWTOA_calc_model(im,:,:)); var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ',im=' num2str(im)];
%                     %         %run plotting script
%                     %         figure
%                     %         UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                     %         lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                     %         caxis([0 250]);
%                     %
%                     %         ioverride_LAT_plots=0;
%                     %         dat_modis = squeeze(SW_up_TOA_dat_ens_monthly(im,:,:)); var_UM = 'Actual SWTOA'; subtitle_str = [var_UM ',im=' num2str(im)];
%                     %         %run plotting script
%                     %         figure
%                     %         UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                     %         lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                     %         caxis([0 250]);
%                     
%                     ioverride_LAT_plots=0;
%                     dat_modis = meanNoNan(surf_albedo(im,:,:),1); var_UM = 'Monthly mean surface albedo'; subtitle_str = [var_UM ' im=' num2str(imon)];
%                     %run plotting script
%                     figure
%                     UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                     lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                     caxis([0 0.2]);
%                     %
%                 end
                
                
%                 %For calculation of map of trends for checking against actual trend map
%                 dat_ukesm.years_ukesm_1d = years_sw_calc;
%                 dat_ukesm.dat_annual = SWTOA_calc_model_annual;                
%                 MIP = 'DAMIP hist-aer SW calculated'; %Be careful with
%                 this as it screws up the MIP value needed in this
%                 calculation script...
%                 var_ukesm = 'SW_up_TOA';
                
                
            end
            
            
            %% Do sensitivity test SW calcs
            
            
            
            %% Keep Cf constant
            
            %istart=1984-1850+1; %starting index to take average over to use as constant value
            
            %Keep cf constant first
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            cf2 = cf;
            
            mean_cf = meanNoNan(cf(inds_trend,:,:),1);
            cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
            cf_const = permute(cf_const,[3 1 2]);
            cf = cf_const;
            
            %keep these as originals
            A_clear = surf_albedo;
            %trans = trans_orig;
            aod = aod_dat_ens_mean;
            aaod = aaod_dat_ens_mean;
            SWclear = sw_TOA_up_clear_dat_ens_mean;
            
            
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %alrady in in-cloud LWP
            end
            
            nd = Nd_dat_ens_mean;
            SWclear = sw_TOA_up_clear_dat_ens_mean;
            
            
            %[Ac_calc_model_cf,tau_calc_model_cf,A_calc_model_cf,SWTOA_calc_model_cf] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
            %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual_cf
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual_cf(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            
            
            %% Keep Nd constant
            
            %keep these as originals for now
            A_clear = surf_albedo;
            %trans = trans_orig;
            aod = aod_dat_ens_mean;
            aaod = aaod_dat_ens_mean;
            SWclear = sw_TOA_up_clear_dat_ens_mean;
            
            
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %alrady in in-cloud LWP
            end
            
            nd = Nd_dat_ens_mean;
            mean_val = meanNoNan(nd(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(nd,1)]);
            val_const = permute(val_const,[3 1 2]);
            nd = val_const;
            
            
            %[Ac_calc_model_Nd,tau_calc_model_Nd,A_calc_model_Nd,SWTOA_calc_model_Nd] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
            %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual_Nd
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual_Nd(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            
            
            %% Keep LWP constant
            
            %keep these as originals for now
            A_clear = surf_albedo;
            %trans = trans_orig;
            aod = aod_dat_ens_mean;
            aaod = aaod_dat_ens_mean;
            SWclear = sw_TOA_up_clear_dat_ens_mean;
                        
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %alrady in in-cloud LWP
            end
            
            mean_val = meanNoNan(lwp(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(lwp,1)]);
            val_const = permute(val_const,[3 1 2]);
            lwp = val_const;
            
            nd = Nd_dat_ens_mean;
            
            %[Ac_calc_model_lwp,tau_calc_model_lwp,A_calc_model_lwp,SWTOA_calc_model_lwp] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);;
            %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual_lwp
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual_lwp(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            
            %% Keep everything constant except LWPic            
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            cf2 = cf;
            
            mean_cf = meanNoNan(cf(inds_trend,:,:),1);
            cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
            cf_const = permute(cf_const,[3 1 2]);
            cf = cf_const;
            
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %already in in-cloud LWP
            end
            
            nd = Nd_dat_ens_mean;
            mean_val = meanNoNan(nd(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(nd,1)]);
            val_const = permute(val_const,[3 1 2]);
            nd = val_const;
            
            mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
            val_const = permute(val_const,[3 1 2]);
            A_clear = val_const;
            
%             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             trans = val_const;
            
            mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aod = val_const;
            
            mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aaod = val_const; 
            
            mean_val = meanNoNan(sw_TOA_up_clear_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(sw_TOA_up_clear_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            SWclear = val_const; 
            
            
            
            
            
            %[Ac_calc_model_lwp_vary,tau_calc_model_lwp_vary,A_calc_model_lwp_vary,SWTOA_calc_model_lwp_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
            %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual_lwp_vary
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual_lwp_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            
            
            %% Keep everything constant except Nd.
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            cf2 = cf;
            
            mean_cf = meanNoNan(cf(inds_trend,:,:),1);
            cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
            cf_const = permute(cf_const,[3 1 2]);
            cf = cf_const; %Prob shouldn't use the constant cf for the LWP calc? Only the case for lwp_type='old', though
            
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %already in in-cloud LWP
            end
            
            mean_val = meanNoNan(lwp(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(lwp,1)]);
            val_const = permute(val_const,[3 1 2]);
            lwp = val_const;
                        
            mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
            val_const = permute(val_const,[3 1 2]);
            A_clear = val_const;
            
%             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             trans = val_const;
            
            mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aod = val_const;
            
            mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aaod = val_const;  
            
            mean_val = meanNoNan(sw_TOA_up_clear_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(sw_TOA_up_clear_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            SWclear = val_const; 
            
            
            nd = Nd_dat_ens_mean; %time-varying value
            %mean_val = meanNoNan(nd(inds_trend,:,:),1);
            %val_const = repmat(mean_val,[1 1 size(nd,1)]);
            %val_const = permute(val_const,[3 1 2]);
            %nd = val_const;
            
            %[Ac_calc_model_cf_lwp,tau_calc_model_cf_lwp,A_calc_model_cf_lwp,SWTOA_calc_model_cf_lwp] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            %[Ac_calc_model_Nd_vary,tau_calc_model_Nd_vary,A_calc_model_cf_lwp,SWTOA_calc_model_Nd_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
            
            clear SWTOA_calc_model_annual_Nd_vary
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                %SWTOA_calc_model_annual_cf_lwp(iy,:,:) = meanNoNan(SWTOA_calc_model_cf_lwp(istart:istart+11,:,:),1);
                SWTOA_calc_model_annual_Nd_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            
            
            
            
            
            %% Keep everything else constant and vary cf.
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            cf2 = cf; %Keep the varying CF for LWPic calc
            
            %mean_cf = meanNoNan(cf(inds_trend,:,:),1);
            %cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
            %cf_const = permute(cf_const,[3 1 2]);
            %cf = cf_const;
            
            switch lwp_type
                case 'old'                    
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %already in in-cloud LWP
            end
            
            mean_val = meanNoNan(lwp(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(lwp,1)]);
            val_const = permute(val_const,[3 1 2]);
            lwp = val_const;
            
            nd = Nd_dat_ens_mean;
            mean_val = meanNoNan(nd(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(nd,1)]);
            val_const = permute(val_const,[3 1 2]);
            nd = val_const;
            
            mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
            val_const = permute(val_const,[3 1 2]);
            A_clear = val_const;
            
%             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             trans = val_const;
            
            mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aod = val_const;
            
            mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aaod = val_const;    
            
            mean_val = meanNoNan(sw_TOA_up_clear_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(sw_TOA_up_clear_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            SWclear = val_const; 
            
            %[Ac_calc_model_cf_vary,tau_calc_model_cf_vary,A_calc_model_cf_vary,SWTOA_calc_model_cf_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
            %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual_cf_vary
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual_cf_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            
%% Keep everything else constant and vary clear_sky flux
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            %cf2 = cf; %Keep the varying CF for LWPic calc
            
            mean_cf = meanNoNan(cf(inds_trend,:,:),1);
            cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
            cf_const = permute(cf_const,[3 1 2]);
            cf = cf_const;
            
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %already in in-cloud LWP
            end
            
            mean_val = meanNoNan(lwp(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(lwp,1)]);
            val_const = permute(val_const,[3 1 2]);
            lwp = val_const;
            
            nd = Nd_dat_ens_mean;
            mean_val = meanNoNan(nd(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(nd,1)]);
            val_const = permute(val_const,[3 1 2]);
            nd = val_const;
            
            mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
            val_const = permute(val_const,[3 1 2]);
            A_clear = val_const;
            %A_clear = surf_albedo; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
            
%             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             trans = val_const;
            
            mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aod = val_const;
            
            mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aaod = val_const;   
            
            SWclear = sw_TOA_up_clear_dat_ens_mean;
            
            %[Ac_calc_model_albedo_vary,tau_calc_model_albedo_vary,A_calc_model_albedo_vary,SWTOA_calc_model_albedo_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
            %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual_clear_sky_vary
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual_clear_sky_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            

            
            
            
%% Keep everything else constant and vary surface albedo
            
            %cf = lowcf_dat_ens_mean;
            cf = totcf_dat_ens_mean;
            %cf2 = cf; %Keep the varying CF for LWPic calc
            
            mean_cf = meanNoNan(cf(inds_trend,:,:),1);
            cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
            cf_const = permute(cf_const,[3 1 2]);
            cf = cf_const;
            
            switch lwp_type
                case 'old'
                    cf2 = cf;
                    cfmin_lwp = cf_min;
                    %cfmin_lwp = 0.05;
                    cf2(cf<cfmin_lwp)=cfmin_lwp;
                    
                    lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
                    
                    
                case 'new'
                    lwp = lwp_dat_ens_mean; %already in in-cloud LWP
            end
            
            mean_val = meanNoNan(lwp(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(lwp,1)]);
            val_const = permute(val_const,[3 1 2]);
            lwp = val_const;
            
            nd = Nd_dat_ens_mean;
            mean_val = meanNoNan(nd(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(nd,1)]);
            val_const = permute(val_const,[3 1 2]);
            nd = val_const;
            
%             mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             A_clear = val_const;
            A_clear = surf_albedo; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
            
%             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             trans = val_const;
            
            mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aod = val_const;
            
            mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            aaod = val_const;    
            
            mean_val = meanNoNan(sw_TOA_up_clear_dat_ens_mean(inds_trend,:,:),1);
            val_const = repmat(mean_val,[1 1 size(sw_TOA_up_clear_dat_ens_mean,1)]);
            val_const = permute(val_const,[3 1 2]);
            SWclear = val_const; 
            
            %[Ac_calc_model_albedo_vary,tau_calc_model_albedo_vary,A_calc_model_albedo_vary,SWTOA_calc_model_albedo_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
            clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
            [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
            %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
            
            clear SWTOA_calc_model_annual_albedo_vary
            for iy=1:length(years_sw_calc)
                istart=(iy-1)*12+1;
                SWTOA_calc_model_annual_albedo_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
            end
            
            
            %% Keep everything else constant and vary transmissivity
%             
%             %cf = lowcf_dat_ens_mean;
%             cf = totcf_dat_ens_mean;
%             %cf2 = cf; %Keep the varying CF for LWPic calc
%             
%             mean_cf = meanNoNan(cf(inds_trend,:,:),1);
%             cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
%             cf_const = permute(cf_const,[3 1 2]);
%             cf = cf_const;
%             
%             switch lwp_type
%                 case 'old'
%                     cf2 = cf;
%                     cfmin_lwp = cf_min;
%                     %cfmin_lwp = 0.05;
%                     cf2(cf<cfmin_lwp)=cfmin_lwp;
%                     
%                     lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%                     
%                     
%                 case 'new'
%                     lwp = lwp_dat_ens_mean; %already in in-cloud LWP
%             end
%             
%             mean_val = meanNoNan(lwp(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(lwp,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             lwp = val_const;
%             
%             nd = Nd_dat_ens_mean;
%             mean_val = meanNoNan(nd(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(nd,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             nd = val_const;
%             
%             mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1); 
%             val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             A_clear = val_const;
% %             A_clear = surf_albedo; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
%             
% %             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
% %             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
% %             val_const = permute(val_const,[3 1 2]);
% %             trans = val_const;
%              trans = trans_orig; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
%             
%             mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             aod = val_const;
%             
%             mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             aaod = val_const; 
%              
%             %[Ac_calc_model_trans_vary,tau_calc_model_trans_vary,A_calc_model_trans_vary,SWTOA_calc_model_trans_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
                %clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
%             [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
%             %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
%             
%             clear SWTOA_calc_model_annual_trans_vary
%             for iy=1:length(years_sw_calc)
%                 istart=(iy-1)*12+1;
%                 SWTOA_calc_model_annual_trans_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
%             end
            
            

%             %% Keep everything else constant and vary AOD
%             
%             %cf = lowcf_dat_ens_mean;
%             cf = totcf_dat_ens_mean;
%             %cf2 = cf; %Keep the varying CF for LWPic calc
%             
%             mean_cf = meanNoNan(cf(inds_trend,:,:),1);
%             cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
%             cf_const = permute(cf_const,[3 1 2]);
%             cf = cf_const;
%             
%             switch lwp_type
%                 case 'old'
%                     cf2 = cf;
%                     cfmin_lwp = cf_min;
%                     %cfmin_lwp = 0.05;
%                     cf2(cf<cfmin_lwp)=cfmin_lwp;
%                     
%                     lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%                     
%                     
%                 case 'new'
%                     lwp = lwp_dat_ens_mean; %already in in-cloud LWP
%             end
%             
%             mean_val = meanNoNan(lwp(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(lwp,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             lwp = val_const;
%             
%             nd = Nd_dat_ens_mean;
%             mean_val = meanNoNan(nd(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(nd,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             nd = val_const;
%             
%             mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1); 
%             val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             A_clear = val_const;
% %             A_clear = surf_albedo; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
%             
% %             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
% %             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
% %             val_const = permute(val_const,[3 1 2]);
% %             trans = val_const;
%              %trans = trans_orig; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
%             
% %             mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
% %             val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
% %             val_const = permute(val_const,[3 1 2]);
% %             aod = val_const;
%             aod = aod_dat_ens_mean;
%             
%             mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             aaod = val_const; 
%              
%             %[Ac_calc_model_trans_vary,tau_calc_model_trans_vary,A_calc_model_trans_vary,SWTOA_calc_model_trans_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
%             clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
%             [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
%             %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
%             
%             clear SWTOA_calc_model_annual_aod_vary
%             for iy=1:length(years_sw_calc)
%                 istart=(iy-1)*12+1;
%                 SWTOA_calc_model_annual_aod_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
%             end
%             
%             
%             %% Keep everything else constant and vary AAOD
%             
%             %cf = lowcf_dat_ens_mean;
%             cf = totcf_dat_ens_mean;
%             %cf2 = cf; %Keep the varying CF for LWPic calc
%             
%             mean_cf = meanNoNan(cf(inds_trend,:,:),1);
%             cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
%             cf_const = permute(cf_const,[3 1 2]);
%             cf = cf_const;
%             
%             switch lwp_type
%                 case 'old'
%                     cf2 = cf;
%                     cfmin_lwp = cf_min;
%                     %cfmin_lwp = 0.05;
%                     cf2(cf<cfmin_lwp)=cfmin_lwp;
%                     
%                     lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%                     
%                     
%                 case 'new'
%                     lwp = lwp_dat_ens_mean; %already in in-cloud LWP
%             end
%             
%             mean_val = meanNoNan(lwp(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(lwp,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             lwp = val_const;
%             
%             nd = Nd_dat_ens_mean;
%             mean_val = meanNoNan(nd(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(nd,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             nd = val_const;
%             
%             mean_val = meanNoNan(surf_albedo(inds_trend,:,:),1); 
%             val_const = repmat(mean_val,[1 1 size(surf_albedo,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             A_clear = val_const;
% %             A_clear = surf_albedo; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
%             
% %             mean_val = meanNoNan(trans_orig(inds_trend,:,:),1);
% %             val_const = repmat(mean_val,[1 1 size(trans_orig,1)]);
% %             val_const = permute(val_const,[3 1 2]);
% %             trans = val_const;
%              %trans = trans_orig; %need to make this the full field and not use inds_trend since this is what is done with cf, etc.
%             
%             mean_val = meanNoNan(aod_dat_ens_mean(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(aod_dat_ens_mean,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             aod = val_const;
% %            aod = aod_dat_ens_mean;
%             
% %             mean_val = meanNoNan(aaod_dat_ens_mean(inds_trend,:,:),1);
% %             val_const = repmat(mean_val,[1 1 size(aaod_dat_ens_mean,1)]);
% %             val_const = permute(val_const,[3 1 2]);
% %             aaod = val_const; 
%             aaod = aaod_dat_ens_mean;
%              
%             %[Ac_calc_model_trans_vary,tau_calc_model_trans_vary,A_calc_model_trans_vary,SWTOA_calc_model_trans_vary] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
%             clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = SWclear;
%             [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
%             %SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
%             
%             clear SWTOA_calc_model_annual_aaod_vary
%             for iy=1:length(years_sw_calc)
%                 istart=(iy-1)*12+1;
%                 SWTOA_calc_model_annual_aaod_vary(iy,:,:) = meanNoNan(SWTOA_calc_model_sens(istart:istart+11,:,:),1);
%             end
%             
                        

            
            
%             % Keep cf and LWPic constant and vary Nd - using CF+LWPic values from specific times (start times of trend period) to test
%             % whether this makes a difference due to cloud fields chagning, etc.
%             
%             cf = lowcf_dat_ens_mean;
%             cf = totcf_dat_ens_mean;
%             cf2 = cf;
%             
%             mean_cf = meanNoNan(cf(inds_trend,:,:),1);
%             mean_cf = meanNoNan(cf(inds_trend(1:12),:,:),1); %use the CF from the start of the trend period
%             cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
%             cf_const = permute(cf_const,[3 1 2]);
%             cf = cf_const;
%             
%             switch lwp_type
%                 case 'old'
%                     cf2 = cf;
%                     cfmin_lwp = cf_min;
%                     cfmin_lwp = 0.05;
%                     cf2(cf<cfmin_lwp)=cfmin_lwp;
%                     
%                     lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%                     
%                     
%                 case 'new'
%                     lwp = lwp_dat_ens_mean; %already in in-cloud LWP
%             end
%             
%             mean_val = meanNoNan(lwp(inds_trend,:,:),1);
%             mean_val = meanNoNan(lwp(inds_trend(1:12),:,:),1); %use the CF from the start of the trend period
%             val_const = repmat(mean_val,[1 1 size(lwp,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             lwp = val_const;
%             
%             nd = Nd_dat_ens_mean;
%             mean_val = meanNoNan(nd(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(nd,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             nd = val_const;
%             
%             [Ac_calc_model_cf_lwp_Pstart,tau_calc_model_cf_lwp_Pstart,A_calc_model_cf_lwp_Pstart,SWTOA_calc_model_cf_lwp_Pstart] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
%             [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
%             SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
%             
%             clear SWTOA_calc_model_annual_cf_lwp_Pstart
%             for iy=1:length(years_sw_calc)
%                 istart=(iy-1)*12+1;
%                 SWTOA_calc_model_annual_cf_lwp_Pstart(iy,:,:) = meanNoNan(SWTOA_calc_model_cf_lwp_Pstart(istart:istart+11,:,:),1);
%             end
% 
%             
% % Keep cf and LWPic constant and vary Nd - using CF+LWPic values from specific times (end time of trend periods)
%             
%             cf = lowcf_dat_ens_mean;
%             cf = totcf_dat_ens_mean;
%             cf2 = cf;
%             
%             mean_cf = meanNoNan(cf(inds_trend,:,:),1);
%             mean_cf = meanNoNan(cf(inds_trend(end-11:end),:,:),1); %use the mean from the last year of the trend period
%             cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
%             cf_const = permute(cf_const,[3 1 2]);
%             cf = cf_const;
%             
%             switch lwp_type
%                 case 'old'
%                     cf2 = cf;
%                     cfmin_lwp = cf_min;
%                     cfmin_lwp = 0.05;
%                     cf2(cf<cfmin_lwp)=cfmin_lwp;
%                     
%                     lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%                     
%                     
%                 case 'new'
%                     lwp = lwp_dat_ens_mean; %already in in-cloud LWP
%             end
%             
%             mean_val = meanNoNan(lwp(inds_trend,:,:),1);
%             mean_val = meanNoNan(lwp(inds_trend(end-11:end),:,:),1); %use the mean from the last year of the trend period
%             val_const = repmat(mean_val,[1 1 size(lwp,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             lwp = val_const;
%             
%             nd = Nd_dat_ens_mean;
%             mean_val = meanNoNan(nd(inds_trend,:,:),1);
%             val_const = repmat(mean_val,[1 1 size(nd,1)]);
%             val_const = permute(val_const,[3 1 2]);
%             nd = val_const;
%             
%             [Ac_calc_model_cf_lwp_Pend,tau_calc_model_cf_lwp_Pend,A_calc_model_cf_lwp_Pend,SWTOA_calc_model_cf_lwp_Pend] = calc_SW2(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu);
%             [Ac_calc_model_sens,tau_calc_model_sens,A_calc_model_sens,SWTOA_calc_model_sens,SWTOA_clear_sky_calc_model_sens] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,cf,lwp*1e3,nd,sw,A_clear,trans_orig);
%             SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);
%             
%             clear SWTOA_calc_model_annual_cf_lwp_Pend
%             for iy=1:length(years_sw_calc)
%                 istart=(iy-1)*12+1;
%                 SWTOA_calc_model_annual_cf_lwp_Pend(iy,:,:) = meanNoNan(SWTOA_calc_model_cf_lwp_Pend(istart:istart+11,:,:),1);
%             end
%             
%             
%                         
%             
            
            %% plot SW
            % var_UM = ['SW TOA model ' str_type_calc ' calc'];
            % dat_modis = SWTOA_calc_model_mean;
            % subtitle_str = var_UM;
            % add_str='';
            % tit_str_clean='Model values';
            % figure
            % ioverride_proj_type=1;
            % proj_type_DRIVER='ortho';
            % irestrict_domain_DRIVER=0;
            % icontour_DRIVER=0;
            % ioverride_LAT_plots=0;
            % ibias_contour=0; cont_col_str_DRIVER='k'; %Colour for the contour lines
            % iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
            % ioverride_ticks_DRIVER=1;
            % icoarse_grain=0;
            % time_round = '';
            % time_format_str=' UTC';
            
            load(load_file,'gcm_Plat2D_UM','gcm_Plon2D_UM','gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM');
            
            irestrict_domain_DRIVER=1;
            ioverride_proj_type=1;
            proj_type_DRIVER='ortho';
            ioverride_LAT_plots=0;
            iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
            ioverride_ticks_DRIVER=1;
            icoarse_grain=0;
            icontour_DRIVER=0; cont_col_str_DRIVER='';
            time_round=0; time_format_str='';
            isave_plot=0;
            iplot_wind_arrows=0;
            
            
            var_UM = ['SW TOA model ' str_type_calc ' calc'];
            dat_modis = SWTOA_calc_model_mean;
            subtitle_str = var_UM;
            
            UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
            lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
            caxis([0 100]);
            xlabel(hc,'W m^{-2}'); %label the colour bar
            
            %Plot the actual SW TOA up
            subtitle_str = ['SW TOA model ACTUAL'];
            dat_modis = meanNoNan(SW_up_TOA_dat_ens_mean,1);
            figure
            ioverride_proj_type=1;
            proj_type_DRIVER='ortho';
            irestrict_domain_DRIVER=0;
            UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
            lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
            caxis([0 100]);
            xlabel(hc,'W m^{-2}'); %label the colour bar
            
            %%
            if isave_plot_SWTOA_calcs_draft==1
                savename=[savedir_date titlenam_driver];
                clear opts
                %        opts.iplot_png=1;
                opts.iplot_eps=1;
                save_map_plot_data; %sets opts values for saving data
                out_file = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
                %        close(gcf);
            end
            
            % %save the map data for later plotting in a subplot plot
            % dat_file = [out_file '.mat'];
            % save(dat_file,'dat_modis','-V7.3');
            
            %% Save the calculated data in a .mat file
            if isingle_ens==1
                ens_str = ['_ens' num2str(iens)];
            else
                ens_str = '_ens_mean_calc_';
            end
            %savedir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';
            savedir = '/home/disk/eos15/d.grosvenor/UM/ACSIS_Nd_trends/';
            
            savename=[savedir 'offline_SW_calcs_time_varying_start_' num2str(yr_start_trend_box) '_end_' num2str(yr_end_trend_box) '_trend_period_' ...
                num2str(iconstant_trend) expt_str ens_str '.mat'];
            savename = remove_character(savename,' ','_');
            
%             vars = {'SW_up_TOA_dat_annual_ens_mean','SWTOA_calc_model_annual','SWTOA_calc_model', ...
%                 'SWTOA_calc_model_annual_Nd','SWTOA_calc_model_annual_cf','SWTOA_calc_model_annual_lwp','SWTOA_calc_model_annual_lwp_vary',...
%                 'SWTOA_calc_model_annual_cf_vary','SWTOA_calc_model_annual_Nd_vary', 'SWTOA_calc_model_annual_cf_lwp_Pstart', 'SWTOA_calc_model_annual_cf_lwp_Pend',...
%                 'SWTOA_calc_model_annual_albedo_vary','SWTOA_calc_model_annual_trans_vary'};
%    

% %vars including the albedo_vary, aod_vary and aaod_vary
%              vars = {'SW_up_TOA_dat_annual_ens_mean','SWTOA_calc_model_annual','SWTOA_calc_model', ...
%                  'SWTOA_clear_sky_calc_model_annual', 'SWTOA_clear_sky_calc_model', ...
%                 'SWTOA_calc_model_annual_Nd','SWTOA_calc_model_annual_cf','SWTOA_calc_model_annual_lwp','SWTOA_calc_model_annual_lwp_vary',...
%                 'SWTOA_calc_model_annual_cf_vary','SWTOA_calc_model_annual_Nd_vary', ...
%                 'SWTOA_calc_model_annual_albedo_vary','SWTOA_calc_model_annual_aod_vary','SWTOA_calc_model_annual_aaod_vary'};           

%vars including the albedo_vary, aod_vary and aaod_vary
             vars = {'SW_up_TOA_dat_annual_ens_mean','SWTOA_calc_model_annual','SWTOA_calc_model', ...
                 'SWTOA_clear_sky_calc_model_annual', 'SWTOA_clear_sky_calc_model', ...
                'SWTOA_calc_model_annual_Nd','SWTOA_calc_model_annual_cf','SWTOA_calc_model_annual_lwp','SWTOA_calc_model_annual_lwp_vary',...
                'SWTOA_calc_model_annual_cf_vary','SWTOA_calc_model_annual_Nd_vary', ...
                'SWTOA_calc_model_annual_clear_sky_vary','SWTOA_calc_model_annual_albedo_vary'};           
                      
            
            for ivar=1:length(vars)
                if ivar==1
                    iappend=0;
                else
                    iappend=1;
                end
                save_vars_mat_func(savename,vars{ivar},eval(vars{ivar}),iappend);
            end
            save(savename,'-APPEND','-V7.3', 'years_sw_calc');
            
            % save(savename,'-V7.3', ...
            % 'years_sw_calc','SW_up_TOA_dat_annual_ens_mean','SW_up_TOA_dat_annual_ens','SWTOA_calc_model_annual','SWTOA_calc_model', ...
            %     'SWTOA_calc_model_annual_Nd','SWTOA_calc_model_annual_cf','SWTOA_calc_model_annual_lwp','SWTOA_calc_model_annual_Nd_cf',...
            %     'SWTOA_calc_model_annual_lwp_Nd','SWTOA_calc_model_annual_cf_lwp', 'SWTOA_calc_model_annual_cf_lwp_Pstart', 'SWTOA_calc_model_annual_cf_lwp_Pend');
            
            
            %zero tot values to start with
            if iens_count==1
                for ivar=1:length(vars)
                    eval_str = [vars{ivar} '_TOT = 0;'];
                    eval(eval_str);
                end
            end
            
            for ivar=1:length(vars)
                eval_str = [vars{ivar} '_TOT = ' vars{ivar} '_TOT + ' vars{ivar} ' ./ Nens;'];
                eval(eval_str);
            end
            %The above will give a NaN value even if just one of the members is NaN,
            %but that should be ok in this case.
            
        end %iens loop
        
        %Now set the TOT array to the ensemble mean and save that too
        if isingle_ens==1
            
            for ivar=1:length(vars)
                eval_str = [vars{ivar} ' = ' vars{ivar} '_TOT;'];
                eval(eval_str);
            end
            
            ens_str = '_post_calc_ens_mean';
            savename=[savedir 'offline_SW_calcs_time_varying_start_' num2str(yr_start_trend_box) '_end_' num2str(yr_end_trend_box) '_trend_period_' ...
                num2str(iconstant_trend) expt_str ens_str '.mat'];
            savename = remove_character(savename,' ','_');
            for ivar=1:length(vars)
                if ivar==1
                    iappend=0;
                else
                    iappend=1;
                end
                save_vars_mat_func(savename,vars{ivar},eval(vars{ivar}),iappend);
            end
            save(savename,'-APPEND','-V7.3', 'years_sw_calc','SW_up_TOA_dat_annual_ens'); %will just save the full ens TOA data
            %in the ens mean file.
            
        end
        
    end %idamip_run loop
    
    
end %period loop

%% Single location calculation for testing.
ido_single_loc=1;
if ido_single_loc==1
    
    lat = 35; lon = -72.5; %off coast of US
    lat = 53; lon = -53; %off coast of US
    
    [minval,ilat] = min(abs(gcm_Plat2D_UM(:,1)-lat));
    [minval,ilon] = min(abs(gcm_Plon2D_UM(1,:)-lon));
    
    surf_albedo_point = surf_albedo(:,ilat,ilon);
    %surf_albedo_point = 0.24;
    
    %Using a constant surf albedo of 0.24 and scaling factor of 0.4333 matches very
    %well!
    
    sw_in_point = sw(:,ilat,ilon);
   % sw_in_point = sw(:,ilat,ilon) ./ me_cos_ALL(:,ilat,ilon); %test using the slantwise component of SWin
    
    
    
    f_upscatter_point = f_upscatter(:,ilat,ilon);
    %f_upscatter_point = 0.7*f_upscatter(:,ilat,ilon);
    %f_upscatter_point = 0.5;
    %f_upscatter_point = 0.25;
    f_upscatter_point = 0.4;
    
    aod_point = aod_dat_ens_mean(:,ilat,ilon);
    %aod_point = (aod_dat_ens_mean(:,ilat,ilon) + dust_aod_dat_ens_mean(:,ilat,ilon));
    %aod_point = aod_point*0.5; %matches very well when using f_up=0.25 and
                %a constant offset of 15 W/m2!
    %aod_point = 0;
    
    aaod_point = aaod_dat_ens_mean(:,ilat,ilon);
    %aaod_point = 20*aaod_dat_ens_mean(:,ilat,ilon);
    %aaod_point = 0;
    
    me_cosSZA_point = me_cos_ALL(:,ilat,ilon);
    me_cosSZA_point = 0.5;
    
    %Using model SW clear-sky :-
    clear opts; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=1; opts.SWclear = sw_TOA_up_clear_dat_ens_mean(:,ilat,ilon);
    %Calculating SW clear-sky and div by cos SZA :-
    %clear opts; opts.idiv_cosSZA=1; opts.cosSZA = me_cosSZA_point; opts.cf_min = cf_min; opts.iprovide_tau=0; opts.i_Liu=1; opts.iprovide_SWclear=0; opts.SWclear = sw_TOA_up_clear_dat_ens_mean(:,ilat,ilon);
    
    %[Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model,SWTOA_clear_sky_calc_model] = calc_SW2(cf,lwp*1e3,nd,sw,surf_albedo,trans_orig,cf_min,0,NaN,i_Liu);
    [Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model,SWTOA_clear_sky_calc_model] = ...
        calc_SW3(opts,aod_point,aaod_point,f_upscatter_point ,...
        trans_clear_sky,cf(:,ilat,ilon),lwp(:,ilat,ilon)*1e3,nd(:,ilat,ilon),sw_in_point,...
        surf_albedo_point,trans_orig);

    SWTOA_calc_model_mean_point = meanNoNan(SWTOA_calc_model,1);
    
    SSA = (aod_point -  aaod_dat_ens_mean(:,ilat,ilon)) ./ aod_point;
    
    clear SWTOA_calc_model_annual_point SWTOA_clear_sky_calc_model_annual_point
    for iy=1:length(years_sw_calc)
        istart=(iy-1)*12+1;
        SWTOA_calc_model_annual_point(iy,:,:) = meanNoNan(SWTOA_calc_model(istart:istart+11,:,:),1);
        SWTOA_clear_sky_calc_model_annual_point(iy,:,:) = meanNoNan(SWTOA_clear_sky_calc_model(istart:istart+11,:,:),1);
        cf_annual(iy,:,:) = meanNoNan(cf(istart:istart+11,:,:),1);
        aod_annual(iy,:,:) = meanNoNan(aod_dat_ens_mean(istart:istart+11,:,:),1);
        aaod_annual(iy,:,:) = meanNoNan(aaod_dat_ens_mean(istart:istart+11,:,:),1);
        fup_annual(iy,:,:) = meanNoNan(f_upscatter(istart:istart+11,:,:),1);
        ssa_annual(iy,:,:) = meanNoNan(SSA(istart:istart+11,:,:),1);
        %dust_aod_annual(iy,:,:) = meanNoNan(dust_aod_dat_ens_mean(istart:istart+11,:,:),1);
        sw_annual(iy,:,:) = meanNoNan(sw(istart:istart+11,:,:),1);
    end
    

fscale = 0.4333;
fscale = 38/27.5;
fscale = 1.1854;
fscale = 1;

clear leg_str
figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,SWTOA_calc_model_annual_point*fscale); hold on; leg_str{1}='Calculated';
plot(years_sw_calc,SW_up_TOA_dat_annual_ens_mean(:,ilat,ilon),'r'); leg_str{2}='Actual';
title('All-sky');
legend(leg_str);

clear leg_str
figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w'); 
plot(years_sw_calc,SWTOA_clear_sky_calc_model_annual_point *fscale); hold on; leg_str{1}='Calculated';
plot(years_sw_calc,sw_TOA_up_clear_dat_annual(:,ilat,ilon),'r'); leg_str{2}='Actual';
title('Clear-sky');
legend(leg_str);

end

%%

iplot_other_vars=0;
if iplot_other_vars==1
    
figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,sw_annual(:,ilat,ilon));
title('SW'); 

figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,cf_annual(:,ilat,ilon));
title('CF');

figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,surf_albedo_annual(:,ilat,ilon));
title('Surface albedo');

figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,aod_annual(:,ilat,ilon));
title('AOD');

figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,dust_aod_annual(:,ilat,ilon));
title('Dust AOD');

figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,aaod_annual(:,ilat,ilon));
title('AAOD');

figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,fup_annual(:,ilat,ilon));
title('fup');

figure; set(gcf,'position',[39 438 1447 420]); set(gcf,'color','w');
plot(years_sw_calc,ssa_annual);
title('SSA');

end



