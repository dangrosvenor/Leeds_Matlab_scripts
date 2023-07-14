%saves the longitude transects calculated by waterVap

saveorload='load';
saveorload='save';

switch saveorload
    case 'load'
        
        model_load={'MODIS','AM3','CAM5'};
        date_range='2007-2010';
        
        season ='DJF';
        season ='ALL';
        
        lon_load='10S';
        lon_load='20S';        
        lon_load='30S';

    case 'save'
        filename_lon_save = ['saved_lon_transects_' var_str '_' time_mean_str '_' gcm_str '_' datestr(now,30)]
        
        lon_saved = xdat;
        lon_transect_saved = ydat;
        description_saved = [titlenam labs(1).l];
        
        
        
        save(filename_lon,'lon_saved','lon_transect_saved','description_saved');
end

