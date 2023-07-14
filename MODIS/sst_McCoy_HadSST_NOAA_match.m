%Convert Dan McCoy's weekly SST dataset to match MODIS arrays
%Is from in-situ and satellite data

sst=load('/home/disk/eos7/dtmccoy/research/SST/SST.mat');

%Matlab time of the dataset
sst.time2 = sst.time+datenum('01-Jan-1800');

[temp,modis_time_matlab] = date_from_day_of_year_func(daynum_timeseries3,modisyear_timeseries3);

ilat = find(sst.lat>lat_restrict(1) & sst.lat <= lat_restrict(2));
lon_restrict2=lon_restrict;
lon_restrict2(lon_restrict<0) = lon_restrict2(lon_restrict<0)+360;
ilon = find(sst.lon>lon_restrict2(1) & sst.lon <= lon_restrict2(2));

sst01 = sst.SST(ilon,ilat,:);   
sst01 = permute(sst01,[3 2 1]);

sst02 = interp1(sst.time2,sst01,modis_time_matlab);