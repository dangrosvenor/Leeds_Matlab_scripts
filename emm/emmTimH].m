fact=1e6*28.97/18;
	logflag=0;
	tit(1).tit='Equilibrium Vapour Mixing Ratio (ppmv)';

    nplots2d=1;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    i2d=3; %tells it are labelling own x axis
    xlabelstr='dN/dt - Number of overshoots per month';
    
    minZ=15.8e3;
    maxZ=17e3;
    ncont=15;
    
    %imaxovr=1;
    maxcovOvr=2.8;
    
    %iminovr=1;
    mincovOvr=2.5;
    
    
    notsame=1;
    
    izovr=1; %flag to set own z axis