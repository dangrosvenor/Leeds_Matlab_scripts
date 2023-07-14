function [lon_new,lat_new]=calc_lat_lon_change_for_dx_dy(lon,lat,dx,dy)
%function [lon_new,lat_new]=calc_lat_lon_change_for_dx_dy(lon,lat,dx,dy)
%lat, lon, dx and dy can be vectors
%dx and dy in metres

%what we really want is to know the path that air would take given an
%initial bearing and windspeed. However, we cannot assume the constant bearing
%or windspeed for long anyway due to likely pressure and Coriolis forces.
%So is proabably best just to use orthogonal lat-lon lines locally to
%calculate a change in lat and lon proportional to the wind speed in each
%direction

%calculate the local scale (conversion between lat lon and km)
%scale distance for 1 degree lat and lon
ddeg=0.1;
%m_lldist gives the distance between sucessive lat lon values in the vector
%so make vectors of e.g. lat(1),lat(1)+ddeg,lat(2),lat(2)+ddeg
%could do a for loop, but might be slower
N=length(lat)*2;
lat2(1:2:N-1)=lat;
lat2(2:2:N)=lat+ddeg;
lon2(1:2:N-1)=lon;
lon2(2:2:N)=lon+ddeg;

x_sc = m_lldist(lon2,[lat lat]);  
y_sc = m_lldist([lon lon],lat2);

%now need to exclude the distances between e.g. lat(1)+ddeg and lat(2)
x_sc2 = 1e3*x_sc(1:2:end)/ddeg;  %and convert from km to m
y_sc2 = 1e3*y_sc(1:2:end)/ddeg;

%so 1 degree equates to x_sc m

lon_new = lon + dx./x_sc2';
lat_new = lat + dy./y_sc2';


