

dx=1e3;
dy=1e3;

clear diff

%for itime=1:size(Times,1)
times=1:11;
times=[1:20];
times=[21:80];
%times=23;

f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

iread_0=0;
iread_1=1;

 for itime2=[1:length(times)]
     itime=times(itime2);
     disp(itime);                 
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    
    lats=[-12 -12.67 -12 -12.67];
    lons=[131 131 131.67 131.67];
    
    [ilats,ilons] = getind_latlon_quick(lat2d.var,lon2d.var,lats,lons,0.1); %the last number should be about one tenth of resolution in km
    
    ilats=[9 19];
    ilons=[156 166];
    
    lat=23;
    lon=149;
    dlat=22;
    dlon=52;
    ilats=[lat-dlat lat+dlat];
    ilons=[lon-dlon lon+dlon];
    
    ilats=[1 201];   %65
    ilons=[1 201];
    
%     for ilat2=1:size(qtotal,2)
%         for ilon2=1:size(qtotal,3)
if itime2==999
    

    if iread_0==1
        qtotal_t0=nc{'QICE'}(itime,:,:,:)+nc{'QSNOW'}(itime,:,:,:)+nc{'QGRAUP'}(itime,:,:,:)+nc{'QVAPOR'}(itime,:,:,:); %total water mixing ratio
        potemp_t0=WRFUserARW(nc,'th',itime);
        z_wrf2_t0=WRFUserARW(nc,'Z',itime); %height profile at position of max total water
        pressure_t0=WRFUserARW(nc,'p',itime); %pressure (mb)
        pressure_t0.var=pressure_t0.var*100; %convert to Pa
        temp_t0=WRFUserARW(nc,'tc',itime); %temperature degC 
    end
    
    mass_t0=0;
    
    for ilat2=min(ilats):max(ilats)
        for ilon2=min(ilons):max(ilons)
            
            i380=find(potemp_t0.var(:,ilat2,ilon2)>380); %indices of all potemps above 380K
            i380=i380(1); %first index
            mass_air = pressure_t0.var(i380:end-1,ilat2,ilon2).*28.97e-3/8.3144./(temp_t0.var(i380:end-1,ilat2,ilon2)+273.15) ...
                .* diff(z_wrf2_t0.var(i380:end,ilat2,ilon2)) *dx*dy;
            mass_t0 = mass_t0 + sum( qtotal_t0(i380:end-1,ilat2,ilon2).*mass_air );                        
            
        end
    end

else
    clear total_water
    
    if iread_1==1
        qvap=nc{'QVAPOR'}(itime,:,:,:); %total water mixing ratio
        qtotal=nc{'QICE'}(itime,:,:,:)+nc{'QSNOW'}(itime,:,:,:)+nc{'QGRAUP'}(itime,:,:,:)+qvap; %total water mixing ratio
        
        potemp=WRFUserARW(nc,'th',itime);
        z_wrf2=WRFUserARW(nc,'Z',itime); %height profile at position of max total water
        pressure=WRFUserARW(nc,'p',itime); %pressure (mb)
        pressure.var=pressure.var*100; %convert to Pa
        temp=WRFUserARW(nc,'tc',itime); %temperature degC 
    end    
    
    for icpt=1:201
        for jcpt=1:201
            [Tmin(icpt,jcpt) imin_cpt(icpt,jcpt)]=min(temp.var(:,icpt,jcpt));
            %pot_cpt(i,j)=potemp.var(imin(i,j),i,j);
            %izmin(i,j)=imin;
        end
    end
    
%    for ilat2=min(ilats):max(ilats)
%        for ilon2=min(ilons):max(ilons)

        ilat2=min(ilats):max(ilats);
        ilon2=min(ilons):max(ilons);

clear deviation dev_big total_water vap_water ice_water_iice total_water_iice vap_water_iice total_water2 vap_water2 total_water3 vap_water3
for iz=75:90   %size(potemp.var,1)-1
            
%            i380=find(potemp.var(iz,min(ilats):max(ilats),min(ilons):max(ilons))>380); %indices of all potemps above 380K
            i380=find(imin_cpt(min(ilats):max(ilats),min(ilons):max(ilons))<iz); %indices of all potemps above 380K
            
            %i380=i380(1); %first index
            %i380=83;
            
            press=pressure.var(iz,min(ilats):max(ilats),min(ilons):max(ilons));
            T=temp.var(iz,min(ilats):max(ilats),min(ilons):max(ilons));
            Z=z_wrf2.var(iz:iz+1,min(ilats):max(ilats),min(ilons):max(ilons));
            qt2=qtotal(iz,min(ilats):max(ilats),min(ilons):max(ilons));
            qv2=qvap(iz,min(ilats):max(ilats),min(ilons):max(ilons));
            
            %iend=size(potemp.var,1);
            %iend=84;
            
            if length(i380)>0
                mass_air = press(i380).*28.97e-3/8.3144./(T(i380)+273.15) ...
                    .* (Z(2,i380) - Z(1,i380))' *dx*dy;
                
                qt=qt2(i380);
                qv=qv2(i380);
                
                val=mean(qt2(:))+3*std(qt2(:));   %threshold for selecting overshoot points - 3 std devs away from mean
                val2=mean(qt2(:))+1*std(qt2(:));   %threshold for selecting icy overshoot points - 1 std devs away from mean
                
                ibig=find(qt>val);                
                iice=find(  qt>val | ( qt-qv > 0.025/f & qt>val2 )  );
                
                ifind=find(f*qt2<prctile(f*qt2(:),40));
                val3=mean(f*qt2(ifind))+5*std(f*qt2(ifind));
                ibig2=find(f*qt>val3);
                
                
%                total_water = qtotal(iz,i380).*mass_air;                        
                total_water(iz,1:length(ibig)) = ( qt(ibig)-mean(qt(:)) ).*mass_air(ibig);  
                vap_water(iz,1:length(ibig))   = ( qv(ibig)-mean(qt(:)) ).*mass_air(ibig);                                
                
                total_water_iice(iz,1:length(iice)) = ( qt(iice)-mean(qt(:)) ).*mass_air(iice);  
                vap_water_iice(iz,1:length(iice))   = ( qv(iice)-mean(qt(:)) ).*mass_air(iice); 
                
                ice_water_iice(iz,1:length(iice))   = ( qt(iice)-qv(iice) ).*mass_air(iice); 
                
                total_water2(iz,1:length(ibig2)) = ( qt(ibig2)-val3/f ).*mass_air(ibig2);  
                vap_water2(iz,1:length(ibig2))   = ( qv(ibig2)-val3/f ).*mass_air(ibig2); 
                
                total_water3(iz,1:length(ibig2)) = ( qt(ibig2)-mean(qt(:)) ).*mass_air(ibig2);  
                vap_water3(iz,1:length(ibig2))   = ( qv(ibig2)-mean(qt(:)) ).*mass_air(ibig2); 

%                med = median(total_water);
%                deviation(iz,1:length(i380)) = total_water - med;
                
%                ibig=find( deviation(iz,:) > maxALL(deviation(iz,:)/4) );
%                dev_big(iz,1:length(ibig)) = deviation(iz,ibig);
%                inc_pos2(itime2) = sum(sum(deviation(iz,ibig)));
                %       end
                %    end
            end
            
end
    
%    inc=total_water-mass_t0
     
      tot_prof = sum(total_water,2);
      inc_tot(itime) = sum(tot_prof);
      
      vap_prof = sum(vap_water,2);
      inc_vap(itime) = sum(vap_prof);
      
      tot_prof_iice = sum(total_water_iice,2);
      inc_tot_iice(itime) = sum(tot_prof_iice);
      
      vap_prof_iice = sum(vap_water_iice,2);
      inc_vap_iice(itime) = sum(vap_prof_iice);
      
      ice_prof = sum(ice_water_iice,2);
      inc_ice(itime) = sum(ice_prof);
      
      
      tot_prof2 = sum(total_water2,2);
      inc_tot2(itime) = sum(tot_prof2);
      
      vap_prof2 = sum(vap_water2,2);
      inc_vap2(itime) = sum(vap_prof2);
      
      
      tot_prof3 = sum(total_water3,2);
      inc_tot3(itime) = sum(tot_prof3);
      
      vap_prof3 = sum(vap_water3,2);
      inc_vap3(itime) = sum(vap_prof3);
      
      
      
    %  ib=find(tot_prof>0.5*std(tot_prof));
    %  inc2(itime2) = sum(tot_prof(ib));

%     dev_prof = sum(deviation,2);
%     dev_big_prof = sum(dev_big,2);
%   
%      inc2(itime2) = sum(dev_prof);
%      inc2_big(itime2) = sum(dev_big_prof);
     
     
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
    
    
         
