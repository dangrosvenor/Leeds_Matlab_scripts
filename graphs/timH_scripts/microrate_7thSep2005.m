%idir=1;
    
	logflag=0;
    fact=1e6*28.97/18;
    
    

    
	
    
    minZ=14.7e3;
    maxZ=22e3;
    
    clines=0; %makes black contour lines appear
    clab=1;
    
    %i2d=2; %tells it to label x axis in km
    
    z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    time=GridDan(idir).t+3;
    
    microcase='vapour source';
    
    switch microcase
    case 'vapour source'
    %pdat(i).p=3600*fact*squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,27:28),3));
    
     pdat(i).p=3600*fact*squeeze( icediagALL(idir).i(izmin:izmax,dumprange,28) );

        
    tit(1).tit='Microphysical Source of Vapour (ppmv/hr)';
    figlab='Microphysical Source of Vap time-height';
    
    iminovr=1;
    mincovOvr=-3;

    imaxovr=1;
    maxcovOvr=3;
    
    sig=2;
    
    case 'source ratio'
    pdat(i).p=300 * squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,27:28),3)) ...
        ./ squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,35:36),3)) ;
    
    tit(1).tit='Microphysical Source Fraction';
    figlab='Microphysical Source Fraction';
    
    iminovr=0;
    mincovOvr=-0.0025;
    mincovOvr=-0.5;

    imaxovr=0;
    maxcovOvr=0.0011;
    maxcovOvr=0.5;
    
	end
    