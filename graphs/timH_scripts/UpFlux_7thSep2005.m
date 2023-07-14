%idir=1;
    
	logflag=0;
    fact=1e6*28.97/18;
    fact=1;
    
    %iminovr=1;
%    mincovOvr=(1e-7);
    mincovOvr=(-7);

    
    %imaxovr=1;
    maxcovOvr=(0.05);
    maxcovOvr=(5E-2);

        
	tit(1).tit='Vapour Flux (ppmv m/s)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	%tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    figlab='Mean Ice time-height';
    
    minZ=14.7e3;
    maxZ=22e3;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    %i2d=2; %tells it to label x axis in km
    
    z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    time=GridDan(idir).t+3;
    
    
  
   pdat(i).p=fact*squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,[25 26]),3));