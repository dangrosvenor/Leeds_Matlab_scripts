%X_flt,Y_flt is the positon based on the WRF d03 domain
%calculated when load the flight data

%now calculate the trajectory of the aircraft - will do this for the 
%grid domain of the WRF grid - then will need to convert to be relative to north

%use these differences like they were u and v velocities for the direction calc
%is distance moved in the time period so is essentially the speed of the aircraft
%the actual time period does not matter
xdiff=[0; diff(X_flt)];
ydiff=[0; diff(Y_flt)];

lat_flt=dat_flt(:,col_lat);
lon_flt=dat_flt(:,col_lon);

%find the lat and lon indices for the WRF domain for the aircraft locations
[ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,lat_flt,lon_flt,0.1);

%will use a function designed for calculating the compass wind direction for WRF
clear dir_flt
for iloc=1:length(ilat)
    dir_flt(iloc)=wind_dir_compass_from_uv_wrf(xdiff(iloc),ydiff(iloc),lat2d,lon2d,ilat(iloc),ilon(iloc),DX,DY);
end

