dlat_start = 10;

minLAT=minALL(lat2d(1).var);
minLAT=floor(minLAT/dlat_start)*dlat_start;
maxLAT=maxALL(lat2d(1).var);
maxLAT=ceil(maxLAT/dlat_start)*dlat_start;
minLON=minALL(lon2d(1).var);
minLON=floor(minLON/dlat_start)*dlat_start;
maxLON=maxALL(lon2d(1).var);
maxLON=ceil(maxLON/dlat_start)*dlat_start;

%dfine_lat=0.2;
%dfine_lon=0.5;

dcoarse_lat=2.5;
dcoarse_lon=5;

dfine_lat=dcoarse_lat/4;
dfine_lon=dcoarse_lon/4;

lats=[minLAT:dcoarse_lat:maxLAT];
lons=[minLON:dcoarse_lon:maxLON];
lats_fine=[floor(minLAT):dfine_lat:ceil(maxLAT)];
lons_fine=[floor(minLON):dfine_lon:ceil(maxLON)];

L_lon = size(lat2d.var,2); %end index of longitude axis
L_lat = size(lat2d.var,1); %end index of longitude axis

hold on
for ilat_draw=1:length(lats)
        [ilat_line ilon_line] = getind_latlon_quick(lat2d.var,lon2d.var,repmat(lats(ilat_draw),[1 length(lons_fine)]),lons_fine,0.2);

	[ilon_line,ilat_line]=add_end_points(ilon_line,ilat_line,L_lon,L_lat);

	ilat_line =(ilat_line - 1)*dy_grid;
	ilon_line =(ilon_line - 1)*dx_grid;	

if length(ilon_line)>1
        ilon_line_fine = ilon_line(1):0.1:ilon_line(end);
        ilat_line_fine = interp1(ilon_line,ilat_line,ilon_line_fine);
else	
	ilon_line_fine = ilon_line;
	ilat_line_fine = ilat_line;
end

	plot(ilon_line_fine,ilat_line_fine,'w-');
end


for ilon_draw=1:length(lons)
        [ilat_line ilon_line] = getind_latlon_quick(lat2d.var,lon2d.var,lats_fine,repmat(lons(ilon_draw),[1 length(lats_fine)]),0.2);

	[ilat_line,ilon_line]=add_end_points(ilat_line,ilon_line,L_lat,L_lon);

	ilat_line =(ilat_line - 1)*dy_grid;
	ilon_line =(ilon_line - 1)*dx_grid;

if length(ilat_line)>1
	ilat_line_fine = ilat_line(1):0.1:ilat_line(end);
	ilon_line_fine = interp1(ilat_line,ilon_line,ilat_line_fine);
else
	ilon_line_fine = ilon_line;
	ilat_line_fine = ilat_line;
end
	
        plot(ilon_line_fine,ilat_line_fine,'w-');
end

disp('Done plotting lat lon lines');




