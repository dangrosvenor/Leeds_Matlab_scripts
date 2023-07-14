lat_roth=-67.57;
lon_roth=-68.13;

[ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,lat_roth,lon_roth,0.1);

domain=2;
domain=3;

switch domain
    case 2
        msize=9;
        plus_size=6;
    case 3
        msize=12;
        plus_size=9;    
        
end

plot(x_grid(ilon),y_grid(ilat),'ko','markerfacecolor','w','markersize',msize);
plot(x_grid(ilon),y_grid(ilat),'r+','markersize',plus_size);


