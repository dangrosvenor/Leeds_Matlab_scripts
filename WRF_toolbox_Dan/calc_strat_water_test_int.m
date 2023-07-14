

dx=1e3;
dy=1e3;

clear diff

%for itime=1:size(Times,1)
times=15:25;
%times=[15 20];

 potemp_grid = [280:1:880];

iread_0=1;
iread_1=1;

 for itime2=[1 4]   %:length(times)   
     itime=times(itime2);
     disp(itime);                 
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    
   % [ilats,ilons] = getind_latlon_quick(lat2d.var,lon2d.var,lats,lons,0.1); %the last number should be about one tenth of resolution in km
    
    ilats=[9 19];
    ilons=[156 166];
    
    lat=23;
    lon=149;
    dlat=20;
    dlon=20;
    ilats=[lat-dlat lat+dlat];
    ilons=[lon-dlon lon+dlon];
    
    ilats=[1 100];
    ilons=[1 201];
    
    ilats=[1 size(p_int_t0(1,:,:),2)];
    ilons=[1 size(p_int_t0(1,:,:),3)];
    
    

%     for ilat2=1:size(qtotal,2)
%         for ilon2=1:size(qtotal,3)
if itime2==1
        mass_t0=0;
    
    for ilat2=min(ilats):max(ilats)
        for ilon2=min(ilons):max(ilons)                        
            
            i380=find(potemp_grid>380); %indices of all potemps above 380K
            i380=i380(1); %first index
            
            inan=isnan(p_int_t0(:,ilat2,ilon2)).*potemp_grid';
            if inan(end)==0
                iend=length(inan);
            else
                iend=find(inan>380);
                iend=iend(1)-1; %the final non-NaN value 
            end
                      
            
            mass_air = p_int_t0(i380:iend-1,ilat2,ilon2).*28.97e-3/8.3144./(tc_int_t0(i380:iend-1,ilat2,ilon2)+273.15) ...
                .* diff(z_int_t0(i380:iend,ilat2,ilon2)) *dx*dy;
            mass_t0 = mass_t0 + sum( q_int_t0(i380:iend-1,ilat2,ilon2).*mass_air );                                        
                                    
            
        end
    end
                
    clear total_water
else
         
    
    for ilat2=min(ilats):max(ilats)
        for ilon2=min(ilons):max(ilons)
            
            i380=find(potemp_grid>380); %indices of all potemps above 380K
            i380=i380(1); %first index
            
            inan=isnan(p_int(:,ilat2,ilon2)).*potemp_grid';
            if inan(end)==0  %if the potemp grid doesn't go beyond the values present
                iend=length(inan);
            else
                iend=find(inan>380);
                iend=iend(1)-1; %the final non-NaN value 
            end
                       
            mass_air = p_int(i380:iend-1,ilat2,ilon2).*28.97e-3/8.3144./(tc_int(i380:iend-1,ilat2,ilon2)+273.15) ...
                .* diff(z_int(i380:iend,ilat2,ilon2)) *dx*dy;
            total_water(ilat2,ilon2) = sum( q_int(i380:iend-1,ilat2,ilon2).*mass_air );                        
            
        end
    end
    
    inc=total_water-mass_t0/(size(q_int,2)*size(q_int,3));
      
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
    
    
         
