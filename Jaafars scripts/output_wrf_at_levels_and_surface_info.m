%Site 1 (Main site used in report as of 13th May). Next to road.
LAT(1)=23.9407;
LON(1)=38.3288;

%Site 2 (beach site)
LAT(2) = 23.9555;
LON(2) = 38.247;

%Site 3 (RASS site)
LAT(3) = 23.9708;
LON(3) = 38.3073;

%Site 4 (test for urban area)
LAT(4)=24.0788;
LON(4)=38.0616;

time=10; %put the correct time index for WRF here (index not GMT)
z0=100; %first height
dz=100; %step in height
ztop=3000; %top height
z_required = [z0:dz:ztop]; %defines the height levels

[ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);



for i=1:length(LAT)
    %BE careful about using the right reference height - the reference here
    %is ABV SEA LEVEL
    Z = WRFUserARW(nca(1).nc,'Z',time,ilat(i),ilon(i));
    terr_level = nc{'HGT'}(:,ilat(i),ilon(i));
    terr_level=terr_level(time);
    Z_T = [terr_level+2 Z]; %add air temp at 2 m
                        
    T = WRFUserARW(nc,'tc',time,ilat(i),ilon(i));
    T = [get_wrf_point_surface(nc,'T2',time,ilat(i),ilon(i))-273.15 T];  
    
    u=WRFUserARW(nc,'u',time,ilat(i),ilon(i));
    v=WRFUserARW(nc,'v',time,ilat(i),ilon(i));
    
    u10 = get_wrf_point_surface(nc,'U10',time,ilat(i),ilon(i));
    v10 = get_wrf_point_surface(nc,'V10',time,ilat(i),ilon(i));
    u=[u10 u];
    v=[v10 v];

    WS = sqrt( u.^2 + v.^2 );
    Z_WS=[terr_level+10 Z];
    
    WD=wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat(i),ilon(i),DX,DY);
    
    fprintf(1,['\n\nSite ' num2str(i) ' Elevation=' num2str(terr_level) 'm LAT=' num2str(LAT(i)) ' LON=' num2str(LON(i))...
        ' i=' num2str(ilon(i)) ' j=' num2str(ilat(i)) ' Time=' Times(time,:)]);
    LU_all=nc{'LU_INDEX'}(time,:);
    LU(i)=LU_all(ilat(i),ilon(i));
    
    TSK_all=nc{'TSK'}(time,:);
    TSK(i)=TSK_all(ilat(i),ilon(i));
    
    IVGTYP_all=nc{'IVGTYP'}(time,:);
    IVGTYP(i)=IVGTYP_all(ilat(i),ilon(i));
    
    ISLTYP_all=nc{'ISLTYP'}(time,:);
    ISLTYP(i)=ISLTYP_all(ilat(i),ilon(i));
    
    VEGFRA_all=nc{'VEGFRA'}(time,:);
    VEGFRA(i)=VEGFRA_all(ilat(i),ilon(i));
    
    TSLB=nc{'TSLB'}(time,:,ilat(i),ilon(i));
    
    fprintf(1,['\nLand Use Index = ' num2str(LU(i))]);
    fprintf(1,['\nSkin temperature = ' num2str(TSK(i)-273.15)]);    
    fprintf(1,['\nDominant vegetation type = ' num2str(IVGTYP(i))]);  
    fprintf(1,['\nDominant soil type = ' num2str(ISLTYP(i))]);      
    fprintf(1,['\nVegetation fraction = ' num2str(VEGFRA(i))]);      
    fprintf(1,['\nSoil temperatures = ' num2str(TSLB-273.15)]);    
    
    fprintf(1,'\nz= %f %f %f %f',0,T(1),WS(1),WD(1));
    
    %interpolate and output results
for iz=1:length(z_required)
    Ti=interp1(Z_T,T,terr_level+z_required(iz));
    WSi=interp1(Z_WS,WS,terr_level+z_required(iz));
    WDi=interp1(Z_WS,WD,terr_level+z_required(iz));    
    fprintf(1,'\nz= %f %f %f %f',z_required(iz),Ti,WSi,WDi);
    
    
end


end
