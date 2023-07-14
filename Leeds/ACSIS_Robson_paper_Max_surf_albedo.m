% Do a calculation of the surface albedo for each ensemble member and all DAMIPs, so can get a map 
% for a mask based on surface albedo (sea-ice).
% Doing is for the annual mean surf albedos since using monthly data
% screens a region that is too large.

MIP = 'DAMIP'; DAMIP_runs = {'DAMIP_hist-aer','DAMIP_hist-GHG','DAMIP_hist-nat'};
isingle_ens=1;

t_inds = [1:165]; 
max_surf = 0;
    
    for idamip_run=1:length(DAMIP_runs)
        
        expt_str = DAMIP_runs{idamip_run};                
        
        for iens=1:4            
            
            var_DAMIP = 'rsuscs';
            switch MIP
                case 'DAMIP'
                    fscale = 1;
                    load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str '_all_' var_DAMIP '.mat'];
                    
                otherwise
                    fscale = 1;
                    model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
            end
            
            mat_obj = matfile(load_file);            
            sw_surf_up_clear_dat_ens_mean = squeeze(fscale *  mat_obj.dat_annual_ens(iens,t_inds,:,:));

            
            
             % --- SW surface downwelling CLEAR-SKY - for calculation of transmissivity and surface albedo ---
            var_DAMIP = 'rsdscs';
            switch MIP
                case 'DAMIP'
                    fscale = 1;
                    load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str '_all_' var_DAMIP '.mat'];
                    
                otherwise
                    fscale = 1;
                    model_str = 'UKESM1'; load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_DAMIP '.mat'];
            end
            
            mat_obj = matfile(load_file);            
            sw_surf_down_clear_dat_ens_mean = squeeze(fscale * mat_obj.dat_annual_ens(iens,t_inds,:,:));
       
            
             surf_albedo  = sw_surf_up_clear_dat_ens_mean ./ sw_surf_down_clear_dat_ens_mean;
             
             max_surf_iens = max(surf_albedo,[],1); %max over all months
             max_surf = max(max_surf,max_surf_iens);
            
        end
        
        
    end
    
    max_surf = squeeze(max_surf);
    
    
    savedir = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/'];
    savename = [savedir 'max_surf_albedo_DAMIP_annual'];
    save(savename,'max_surf','-V7.3');
    
    