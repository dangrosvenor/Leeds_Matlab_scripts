function dir=wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat,ilon,DX,DY)
%latest wind direction calcualtor as of 17th Feb 2010
%function dir=wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat,ilon,DX,DY)

xinds = 1:size(lat2d(1).var,2); %lat2d is arranged in LAT,LON order so x is second index (LON)
yinds = 1:size(lat2d(1).var,1);

dx_grid = DX;
dy_grid = DY;

grid_y = (yinds - 1)*dy_grid;
grid_x = (xinds - 1)*dx_grid;

sp=1;

if ilat+sp>size(lat2d.var,1)
    jnorth = ilat - sp;
else
    jnorth = ilat + sp;  %get a reference point to figure out the local direction of the lat/lon lines
end

lons_north = lon2d.var(jnorth,:);

%%%OR 1D interpolate the x and y position along the line

xref = interp1(lons_north,grid_x,lon2d.var(ilat,ilon));
yref = grid_y(jnorth);

%[temp inorth] = min( abs(lons_north - lon2d.var(ilat,ilon) ) );

%or could do contour to find the angle of the same latitude contour

%angle of the local north line relative to the grid
%%%thetaN = atan ( (inorth - ilon) / (jnorth - ilat) );

thetaN = atan ( (xref - grid_x(ilon)) / (yref - grid_y(ilat)) );



dir=wind_dir_from_uv(u,v);

dir = dir - thetaN*180/pi; 
% take away thetaN to give direction relative to north
i360 = find(dir>=360);
dir(i360) = dir(i360) - 360;

i0 = find(dir<0);
dir(i0) = dir(i0) + 360;



