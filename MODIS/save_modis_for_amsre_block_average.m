modis_for_amsre_filename = '/home/disk/eos5/d.grosvenor/AMSRE/modis_saved_data_for_amsre_block_average.mat';
LAT_val = lat_restrict;
LON_val = lon_restrict;
save(modis_for_amsre_filename,'modisyear_timeseries3','daynum_timeseries3','Cloud_Fraction_Liquid','LAT_val','LON_val','years_amsre_str');