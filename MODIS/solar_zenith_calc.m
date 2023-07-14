function solar_zenith = solar_zenith_calc(lat,lon,scantime_matlab)


for i1=1:length(lat(:))
   time_mod = datestr(scantime_matlab(i1));
   location.latitude = lat(i1);
   location.longitude = lon(i1);
   location.altitude = 0;
   sun = sun_position(time_mod,location);
   solar_zenith(i1) = sun.zenith;
end

solar_zenith = reshape(solar_zenith,size(lat));
   
   