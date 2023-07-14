    
	fact=1;
	logflag=0;
    clab=1;
    figlab='max sat ice total';
    izovr=1;
	
    
    %iminovr=[1];
    mincovOvr=[0];
    i2d=0;
    
    
    minZ=15e3;
    maxZ=22e3;
    
    [izmin izmax]=findheight(zmpc,minZ,maxZ);
    
    ncont=40;

    tit(1).tit='Temp (^oC)';
    %tit(1).tit=['Max Vapour in Updraughts (g/kg)'];
    %tit(1).tit=['Max W (m/s)'];
    %tit(1).tit=['Min Temp (^oC)'];
	%tit(2).tit=['Mean Ice Mixing Ratio Time Height for MPC (g/kg)'];
    
    
    timesTH(i).t=TimeMPC(1:end)';
    zz(i).z=zmpc(izmin:izmax);
        
    pdat(i).p=squeeze(min(tmpmpc(:,:,izmin:izmax),[],2))'-273.15;
    
