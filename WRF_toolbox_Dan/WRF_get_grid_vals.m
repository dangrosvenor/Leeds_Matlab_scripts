function [pdat_out,rho_out,lat,lon,wind,ix,iy,ix2,iy2] = WRF_get_grid_vals(x,y,x_grid,y_grid,pdat_orig,rho,lat2d,lon2d,x_grid_orig,zz,u_quiver,v_quiver,ih_wrf)

%displays a few useful things for a given x,y point on a WRF grid performed
%using PlotTime with the 'wrf_plot' plotting case.

fprintf(1,'\nx = %f km',x);
fprintf(1,'\ny = %f km',y);

ix=findheight_nearest(x_grid,x);
iy=findheight_nearest(y_grid,y);

pdat_out = pdat_orig(iy,ix);
fprintf(1,'\npdat = %f ',pdat_out);

rho_out = rho(iy,ix);
fprintf(1,'\nrho on level %d = %f kg m^{-3}',ih_wrf,rho_out);

lat = lat2d.var(iy,ix);
lon = lon2d.var(iy,ix);

fprintf(1,'\nlat = %f deg',lat);
fprintf(1,'\nlon = %f deg',lon);


%now for variables stored on the zoomed in grid.
ix2=findheight_nearest(x_grid_orig,x);
iy2=findheight_nearest(zz(1).z,y);

wind=sqrt(u_quiver(iy2,ix2).^2 + v_quiver(iy2,ix2).^2);

fprintf(1,'\nWind speed on level %d = %f m s^{-1}',ih_wrf,wind);

fprintf(1,'\n');




