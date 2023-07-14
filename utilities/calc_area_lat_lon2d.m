function [area] = calc_area_lat_lon2d(lat2d,lon2d)
%order of lat2d and lon2d should be (ilat,ilon)
%lon values should also start at 0 and go eastwards

%convert to 0 to 360 style
i0 = find(lon2d<0);
lon2d(i0) = lon2d(i0) + 360;

for ilat=1:size(lat2d,1)-1
    for ilon=1:size(lon2d,2)-1
        
        dlat = lat2d(ilat+1,ilon) - lat2d(ilat,ilon);
        dlon = lon2d(ilat,ilon+1) - lon2d(ilat,ilon);
        %mid points - use these to avoid issues at lat=90
        lat = 0.5 * (lat2d(ilat+1,ilon) + lat2d(ilat,ilon));
        lon = 0.5 * (lon2d(ilat,ilon+1) + lon2d(ilat,ilon));
        
   
        dist_lat = distlatlon(lat-dlat/2,lon,lat+dlat/2,lon);
        dist_lon = distlatlon(lat,lon-dlon/2,lat,lon+dlon/2);
        area(ilat,ilon) = dist_lat .* dist_lon;                        
        
    end
end