function [lat,lon]=get_lat_lon_wrf_from_xy(lat2d,lon2d,X,Y,x,y)

ix=findheight_nearest(X,x);
iy=findheight_nearest(Y,y);

lat=lat2d.var(iy,ix);
lon=lon2d.var(iy,ix);