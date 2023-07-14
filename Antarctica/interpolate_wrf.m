LAT=[-67.5702 -67.1420];
LON=[-68.1297 -61.6650];
[ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);

time=16; %44=19, 45=20 UTC 16=21 UTC for ncep

heights=[27.5 78.6];

for i=1:length(ilat)
	z_wrf = WRFUserARW(nc,'Z',time,ilat(i),ilon(i));
	p_wrf = WRFUserARW(nc,'p',time,ilat(i),ilon(i));
	T_wrf = WRFUserARW(nc,'tc',time,ilat(i),ilon(i));
	qv_wrf = nc{'QVAPOR'}(time,:,ilat(i),ilon(i));

	u=WRFUserARW(nc,'u',time,ilat(i),ilon(i));
        v=WRFUserARW(nc,'v',time,ilat(i),ilon(i));        
        wind_wrf = sqrt( u.^2 + v.^2 );

 	u=WRFUserARW(nc,'umeta',time,ilat(i),ilon(i));
        v=WRFUserARW(nc,'vmeta',time,ilat(i),ilon(i));     
        [th,r]=cart2pol(u,v); %th : minus indicates angles anti-clockwise from 0 to 180

        wdir_wrf = (360 + th*180/pi);

	p2(i) = interp1(z_wrf,p_wrf,heights(i),[],'extrap');
	T2(i) = interp1(z_wrf,T_wrf,heights(i),[],'extrap');
	qv2_wrf2(i) = interp1(z_wrf,qv_wrf,heights(i),[],'extrap');	
	wind2(i) = interp1(z_wrf,wind_wrf,heights(i),[],'extrap');	
	wdir2(i) = interp1(z_wrf,wdir_wrf,heights(i),[],'extrap');	


end

%Tsurf=WRFUserARW(nc,'tc',time,1,:,:);

