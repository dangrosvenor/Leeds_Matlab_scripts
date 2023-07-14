%what we really want is to know the path that air would take given an
%initial bearing and windspeed. However, we can assume the constant bearing
%or windspeed for long anyway due to likely pressure and Coriolis forces.
%So is proabably best just to use orthogonal lat-lon lines locally to
%calculate a change in lat and lon proportional to the wind speed in each
%direction

sat_time = 21*3600+15*60; %satellite image time in seconds

dt = sat_time - mpace_time*3600;
%so if the flight time is after the satellite time we want negative dt
%(back trajectory)

windspeed=11.8; %m/s
winddir=67.25; %degrees from north

u = windspeed * sin(winddir*pi/180);
v = windspeed * cos(winddir*pi/180);

[mpace_lon_mapped,mpace_lat_mapped]=calc_lat_lon_change_for_dx_dy(mpace_lon,mpace_lat,u*dt,v*dt);

