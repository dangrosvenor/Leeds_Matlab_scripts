function dat=L2L3_concatenate_daily_arrays(load_filename,var_str,lat_lon_inds_str)
%load the data
%lat_lon_inds_str is a string - e.g. '(1:10,20:30,:)' to tell it what
%indices to output
eval(['load(load_filename,''' var_str ''');']);

dat = eval([var_str '.timeseries3' lat_lon_inds_str ';']);