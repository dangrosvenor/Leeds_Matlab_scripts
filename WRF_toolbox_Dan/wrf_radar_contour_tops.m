

dx=1e3;
dy=1e3;

clear diff

%for itime=1:size(Times,1)
times=218:234;
times=173:189;
times=278:294;
times=270:273;
times=[161:172];
times=[421:500];
%times=226;


%clear maxh iabove

for itime2=1:length(times)   
     itime=times(itime2);
     disp(itime);                        

    conts=[10 20 35 40];
    
    izs=[76:90];
    clear z_slice pot_slice
    for iz=izs
        z_slice(iz,:,:)=WRFUserARW(nc,'Z',time,iz);
        pot_slice(iz,:,:)=WRFUserARW(nc,'th',time,iz);
    end
           
    
    clear z380
    nh=length(z_slice(1,:));
    for ih=1:nh;
        i380=find(pot_slice(:,ih)>380);
        i380=i380(1);
        z380(ih)=z_slice(i380,ih);
    end
    
    z380=median(z380); %take the most common height for 380K to represent the NCEP mesoscale 380K level
    z380_itime(itime)=z380;
    
    clear maxh_slice   
    for iz=izs
        Z=WRFRadarRefl_smallMEM_2(nc,itime,'thompson',iz); %reflectivity for this level
        Z=real(10.*log10(Z)); 
        Z(find(Z(:)<0))=NaN;
        
      for iconts=1:length(conts)        
        cont=conts(iconts);         
        
        idbz=find(Z>cont);
        if length(idbz)>0
            maxh_slice(iz,iconts)=max(z_slice(iz,idbz)); %max height of all the >cont dbz points
        else
            maxh_slice(iz,iconts)=0;
        end               
        
      end
      
      
    end
    
   for iconts=1:length(conts)   
        maxh(itime,iconts)=max( maxh_slice(:,iconts) );
        
        istrat=find(maxh(itime,iconts)>z380);
        if length(istrat)>0
            iabove(itime,iconts)=1;
        else
            iabove(itime,iconts)=0;
        end
        
   end
    
                        
    

end

disp('Finished');
    
    
         
