%% Also calculate the mean values over the different regions.
ilat=0;
for iregion=1:length(iregions)    
    ireg=find(i2==iregion);
    dat_reg=[];
    dat_std_reg=[];
    dat_N_reg=[];
    dat_me_reg=[];
    
    for i=1:length(ireg)                
        ilat = ilat + 1;          
        LAT_val_DRIVER = LATs{ilat}; LON_val_DRIVER = LONs{ilat};
    
        ilat2 = find(MLAT >= LAT_val_DRIVER(1) & MLAT <= LAT_val_DRIVER(2));
        ilon2 = find(MLON >= LON_val_DRIVER(1) & MLON <= LON_val_DRIVER(2));  
        
        dat = dat_modis(ilat2,ilon2);        
        dat_reg = cat(1,dat_reg,dat(:));
        if exist('stdev')
            %stdev, me and nnums calculated in calling script
            dat_std = stdev(ilat2,ilon2);
            dat_std_reg = cat(1,dat_std_reg,dat_std(:));
            dat_N = nnums(ilat2,ilon2);
            dat_N_reg = cat(1,dat_N_reg,dat_N(:));            
            dat_me = me(ilat2,ilon2);
            dat_me_reg = cat(1,dat_me_reg,dat_me(:));                        
        else
            dat_std_reg = NaN;
            dat_N_reg = NaN;
            dat_me_reg = NaN;
        end
        
    end
    
    mean_reg_val(iregion) = meanNoNan(dat_reg,1); %mean of dat_modis data - for std plots this is
    %the mean over all locations of the relative std devs of errors (taken over time)
    max_reg_val(iregion) = maxALL(dat_reg);
    min_reg_val(iregion) = minALL(dat_reg);    
    
    mean_abs_err_reg_val(iregion) = meanNoNan(dat_me_reg,1); %whereas this is the mean of what 
    % was used for the me caluculation (i.e. either relative or abs error)
    mean_std_abs_err_reg_val(iregion) = meanNoNan(dat_std_reg,1); %as above, but for the std dev of
    % errors - i.e. spatial mean of std dev (over time) of errors
    
    %The above count each region equally no matter how many datapoints go
    %into the timeseries
    
    %This weights each regions result by the number of points in their timeseries
    [std_reg_val(iregion), me_reg_val(iregion)]=std_combine2(dat_std_reg',dat_me_reg',dat_N_reg');
end







% General instructions for case 'Generic plot specified outside of script'
%Run plot_global_maps_defaults (might want to change these)
            %Then need to set e.g. P, gcm_Plat2D_CERES, gcm_Plon2D_CERES,
            % gcm_Plat2D_edges_CERES, gcm_Plon2D_edges_CERES],
            % gcm_str, daynum_timeseries3_CERES = 1; gcm_time_UTC_CERES = 1; month_amsre=1; year_amsre=1;
            % mod_data_type='AMSRE';
            % Can set these times to one, or do properly depending on what
            % time screening you want, etc.