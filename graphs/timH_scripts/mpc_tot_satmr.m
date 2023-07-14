    
	fact=1;
	logflag=0;
    clab=1;
    figlab='max sat ice total';
    izovr=1;
	
    
    iminovr=[1];
    mincovOvr=[0];
    i2d=0;
    
    
    minZ=15e3;
    maxZ=22e3;
    
    [izmin izmax]=findheight(zmpc,minZ,maxZ);
    
    ncont=40;

    tit(1).tit=['Max Total Water Sat Ratio wrt ice'];    
    %tit(1).tit=['Max Vapour in Updraughts (g/kg)'];
    %tit(1).tit=['Max W (m/s)'];
    %tit(1).tit=['Min Temp (^oC)'];
	%tit(2).tit=['Mean Ice Mixing Ratio Time Height for MPC (g/kg)'];
    
    
    timesTH(i).t=TimeMPC(1:end)';
    zz(i).z=zmpc(izmin:izmax);

       % icesat=SatVapPress(tmpmpc(:,:,izmin:izmax),'goff','ice',100*prsmpc(:,:,izmin:izmax),1); %ppmv
           %rh=wmvmpc(:,:,izmin:izmax)./icesat;

    %icesat=SatVapPress(tmpmpc(:,:,izmin:izmax),'goff','ice'); %Pa
    %p_i=100*prsmpc(:,:,izmin:izmax).*wmvmpc(:,:,izmin:izmax)./(0.622+wmvmpc(:,:,izmin:izmax));
    %rh=p_i./icesat;
        
    pdat(i).p=squeeze(max(rh,[],2))';
    
    nsig=1;