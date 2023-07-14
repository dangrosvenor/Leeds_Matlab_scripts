

dx=1e3;
dy=1e3;

clear diff

%for itime=1:size(Times,1)
times=15:25;
%times=[15 20];

 for itime2=[1 4]   %:length(times)   
     itime=times(itime2);
     disp(itime);
    
    qtotal=nc{'QICE'}(itime,:,:,:)+nc{'QSNOW'}(itime,:,:,:)+nc{'QGRAUP'}(itime,:,:,:)+nc{'QVAPOR'}(itime,:,:,:); %total water mixing ratio
    potemp=WRFUserARW(nc,'th',itime);
    z_wrf2=WRFUserARW(nc,'Z',itime); %height profile at position of max total water
    pressure=WRFUserARW(nc,'p',itime); %pressure (mb)
    pressure.var=pressure.var*100; %convert to Pa
    temp=WRFUserARW(nc,'tc',itime); %temperature degC 
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    
    [ilats,ilons] = getind_latlon_quick(lat2d.var,lon2d.var,lats,lons,0.1); %the last number should be about one tenth of resolution in km
    
    ilats=[9 19];
    ilons=[156 166];
    
    lat=23;
    lon=149;
    ilats=[lat-5 lat+5];
    ilons=[lon-5 lon+5];
    
%    ilats=[101 201];
%    ilons=[101 201];
    
    total_water=0;
%     for ilat2=1:size(qtotal,2)
%         for ilon2=1:size(qtotal,3)

     for ilat2=min(ilats):max(ilats)
         for ilon2=1:min(ilons):max(ilons)
            
            i380=find(potemp.var(:,ilat2,ilon2)>380); %indices of all potemps above 380K
            i380=i380(1); %first index
            mass_air = pressure.var(i380:end-1,ilat2,ilon2).*28.97e-3/8.3144./(temp.var(i380:end-1,ilat2,ilon2)+273.15) ...
                .* diff(z_wrf2.var(i380:end,ilat2,ilon2)) *dx*dy;
            total_water = total_water + sum( qtotal(i380:end-1,ilat2,ilon2).*mass_air );                        
            
        end
    end
    
    if itime2==1
        tot_inc(itime2) = 0;               
        mass_t0 = total_water;
    else
        tot_inc(itime2) = total_water - mass_t0;
        mass_tot(itime2) = total_water;
    end
    
%     %quick estimate for comparison
%     i380=find(potemp.var(:,ilat(1),ilon(1))>380);
%     mass_air = pressure.var(i380:end-1,ilat(1),ilon(1)).*28.97e-3/8.3144./(temp.var(i380:end-1,ilat(1),ilon(1))+273.15) ...
%                 .* diff(z_wrf2.var(i380:end,ilat(1),ilon(1))) *dx*dy;
%             
%             total_water=0;
%             for ilat2=1:size(qtotal,2)
%                 for ilat2=1:size(qtotal,3)
%                                        
%                     total_water = total_water + sum( qtotal(i380:end-1,ilat2,ilon2).*mass_air );
%                     
%                     if itime==1
%                         tot_inc2(itime) = 0;
%                         mass_t02 = total_water;
%                     else
%                         tot_inc2(itime) = total_water - mass_t0;
%                     end
%                     
%                 end
%             end

end

disp('Finished');
    
    
         
