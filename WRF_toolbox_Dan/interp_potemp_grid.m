iload_t0=0;

if iload_t0==1
    itime=15;
    
    qtotal=nc{'QICE'}(itime,:,:,:)+nc{'QSNOW'}(itime,:,:,:)+nc{'QGRAUP'}(itime,:,:,:)+nc{'QVAPOR'}(itime,:,:,:); %total water mixing ratio
    potemp=WRFUserARW(nc,'th',itime);    
    z_wrf2=WRFUserARW(nc,'Z',itime); %height profile at position of max total water
    pressure=WRFUserARW(nc,'p',itime); %pressure (mb)
    pressure.var=pressure.var*100; %convert to Pa
    temp=WRFUserARW(nc,'tc',itime); %temperature degC 
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    [ilats,ilons] = getind_latlon_quick(lat2d.var,lon2d.var,lats,lons,0.1); %the last number should be about one tenth of resolution in km
    
    total_water=0;
    %     for ilat2=1:size(qtotal,2)
    %         for ilon2=1:size(qtotal,3)
    potemp_grid = [minALL(potemp.var):1:maxALL(potemp.var)];
    
    ilats=[1 101];
    ilons=[101 201];
    
    clear q_int_t0 th_int_t0 tc_int_t0 p_int_t0 z_int_t0;
    for ilat2=min(ilats):max(ilats)
        for ilon2=min(ilons):max(ilons)
            [B,I]=unique(potemp.var(:,ilat2,ilon2));
            
            th_int_t0( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = potemp_grid;
            q_int_t0( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,qtotal(I,ilat2,ilon2),potemp_grid);            
            tc_int_t0( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,temp.var(I,ilat2,ilon2),potemp_grid);
            p_int_t0( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,pressure.var(I,ilat2,ilon2),potemp_grid);
            z_int_t0( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,z_wrf2.var(I,ilat2,ilon2),potemp_grid);            
            
        end
    end
    
end
    

itime=18;
    
    qtotal=nc{'QICE'}(itime,:,:,:)+nc{'QSNOW'}(itime,:,:,:)+nc{'QGRAUP'}(itime,:,:,:)+nc{'QVAPOR'}(itime,:,:,:); %total water mixing ratio
    potemp=WRFUserARW(nc,'th',itime);    
    z_wrf2=WRFUserARW(nc,'Z',itime); %height profile at position of max total water
	pressure=WRFUserARW(nc,'p',itime); %pressure (mb)
	pressure.var=pressure.var*100; %convert to Pa
	temp=WRFUserARW(nc,'tc',itime); %temperature degC 
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    [ilats,ilons] = getind_latlon_quick(lat2d.var,lon2d.var,lats,lons,0.1); %the last number should be about one tenth of resolution in km
    
    total_water=0;
    %     for ilat2=1:size(qtotal,2)
    %         for ilon2=1:size(qtotal,3)
    potemp_grid = [minALL(potemp.var):1:maxALL(potemp.var)];
    
    ilats=[1 101];
    ilons=[101 201];
    
	clear q_int th_int T_int p_int z_int;
    for ilat2=min(ilats):max(ilats)
        for ilon2=min(ilons):max(ilons)
            [B,I]=unique(potemp.var(:,ilat2,ilon2));
            
            q_int( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,qtotal(I,ilat2,ilon2),potemp_grid);
            th_int( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,potemp.var(I,ilat2,ilon2),potemp_grid);
            tc_int( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = 273.15 + interp1(B,temp.var(I,ilat2,ilon2),potemp_grid);
            p_int( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,pressure.var(I,ilat2,ilon2),potemp_grid);
            z_int( :,ilat2-min(ilats)+1,ilon2-min(ilons)+1 ) = interp1(B,z_wrf2.var(I,ilat2,ilon2),potemp_grid);            
            
        end
    end
    
    disp('done potemp interpolation');