function inds=get_inds_constant_lat(lat,lat2d,lon2d,lon_min,lon_max)
% function inds=get_inds_constant_lat(lat,lat2d,lon2d,lon_min,lon_max)
% Can choose the lenght of the slice needed in terms of lon (from lon_min to lon_max)
% for all lons just leave out lon_min and lon_max
% for a constant lon line just swap around input of lat and lon
% perhaps should do it as a range of latitudes - e.g. all points from latA to latB to get a 'bin' of latitude

dfli = 0.05;
if nargin==5
    inds = find(lat2d.var>lat-dfli & lat2d.var<lat+dfli & lon2d.var>lon_min & lon2d.var<lon_max);
else
    inds = find(lat2d.var>lat-dfli & lat2d.var<lat+dfli);
end