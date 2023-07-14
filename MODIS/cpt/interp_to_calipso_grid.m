%CALIPSO data is 90x180 - i.e. 2x2 deg. So interpolate the
%MODIS data to the same grid

MLON_grid = eval(['[floor(minALL(gcm_Plon2D_CALIPSO_monthly)) : 2 : ceil(maxALL(gcm_Plon2D_CALIPSO_monthly)) ];']);
MLAT_grid = eval(['[floor(minALL(gcm_Plat2D_CALIPSO_monthly)) : 2 : ceil(maxALL(gcm_Plat2D_CALIPSO_monthly)) ];']);
[X,Y]= meshgrid(MLON_grid,MLAT_grid);

dlon = mean(diff(MLON_grid)); %this will give the order of the array too, so could be negative
dlat = mean(diff(MLAT_grid)); %if runs from e.g. 90 to -90

MLON_grid_edges = MLON_grid(1) - dlon/2 : dlon : MLON_grid(end) + dlon/2;
MLAT_grid_edges = MLAT_grid(1) - dlat/2 : dlat : MLAT_grid(end) + dlat/2;
[X2,Y2]= meshgrid(MLON_grid_edges,MLAT_grid_edges);

%  Then execute something like this to do the actual interpolation:-
%            modis_cf_grid = griddata(MLON,MLAT,P2,X,Y);  
% where MLAT and MLON are the grid for the data to be interpolated (P2)