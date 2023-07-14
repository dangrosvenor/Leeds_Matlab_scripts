qf=nc{'QICE'}(time,:,:,:)+nc{'QVAPOR'}(time,:,:,:);
[x,I]=maxALL(f*qf(iz,:,:));


clear trop

multi=0;

if multi==0
    z_wrf2=WRFUserARW(nc,'Z',time,I(1),I(2)); %height profile at position of max total water
    temp=WRFUserARW(nc,'tc',time,I(1),I(2));
    pressure=WRFUserARW(nc,'p',time,I(1),I(2)); %pressure in mb
    pressure.var=pressure.var*100;
    %[trop,trop_press,trop2,trop_press2]=wmo_trop(H,temp2,press,smooth,windowSize)
    [trop,trop_press,trop2,trop_press2]=wmo_trop(z_wrf2,temp,pressure.var,0,0);
    
    potemp=WRFUserARW(nc,'th',time);
    i380=find(potemp.var(:,I(1),I(2))>380); %indices of all potemps above 380K
    trop380=z_wrf2(i380(1)); %first height above 380 K
else
    
    for irand=1:10
        
        I(1)=round(rand(1)*size(qf,2));
        I(2)=round(rand(1)*size(qf,3));
        
        
        %for ilat2=1:size(qf,2)
        %    for ilon2=1:size(qf,2)
        %        I(1)=ilat2;
        %        I(2)=ilon2;
        z_wrf2=WRFUserARW(nc,'Z',time,I(1),I(2)); %height profile at position of max total water
        temp=WRFUserARW(nc,'tc',time,I(1),I(2));
        pressure=WRFUserARW(nc,'p',time,I(1),I(2)); %pressure in mb
        pressure.var=pressure.var*100;
        %[trop,trop_press,trop2,trop_press2]=wmo_trop(H,temp2,press,smooth,windowSize)
        [trop(irand),trop_press,trop2,trop_press2]=wmo_trop(z_wrf2,temp,pressure.var,0,0);
        
        %    end
        %end
        
        
    end
    
    
end

        
        