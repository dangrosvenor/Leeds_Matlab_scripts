%works out x and y as a function of a regular lat lon grid

%run latlon_make_fine_grid.m first to make the fine lat lon grid

max_lat=maxALL(lat2d_fine);
min_lat=minALL(lat2d_fine);
dlat=( max_lat-min_lat )/100;
lats=[min_lat:dlat:max_lat];

max_lon=maxALL(lon2d_fine);
min_lon=minALL(lon2d_fine);
dlon=( max_lon-min_lon )/100;
lons=[min_lon:dlon:max_lon];

[LATS,LONS]=meshgrid(lats,lons); %makes lat and lon matrices so that LATS(:),LONS(:) would
%be all the lat,lon pairings

LATS2=LATS(:); %make into a vector
LONS2=LONS(:); %make into a vector

[ilat,ilon] = getind_latlon_quick(lat2d_fine,lon2d_fine,LATS2,LONS2,0.1);

dx = x_fine(2)-x_fine(1);
X_latlon = (reshape(ilon,size(LONS))-1) * dx;

dy = y_fine(2)-y_fine(1);
Y_latlon = (reshape(ilat,size(LATS))-1) * dy;

%filename='.....';
%save(filename,'X_latlon','Y_latlon','LATS','LONS')