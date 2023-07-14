fact=1e6*28.97/18;
	logflag=0;
	tit(1).tit='Vapour Mixing Ratio (ppmv)';
		%tit(1).tit='Vapour Mixing Ratio (ppmv)';
    tit(2).tit='Potential Temperature (K)';
    nplots2d=2;
    
    clines=1; %makes black contour lines appear
    clab=1;
    
    i2d=2; %tells it to label x axis in km

    ncont=15;
    
    %imaxovr=1;
    maxcovOvr=2.8;
    
    %iminovr=1;
    mincovOvr=2.5;
    
    z=GridDan(1).Z;
    
    notsame=1;
    
    itimestamp=1;
    

	switch i
            
        case 1
        [izmin izmax]=findheight(Grid.Z,14e3,22e3);
        timesTH(1).t=Grid.Y1'/1000;
        zz(1).z=GridDan(1).Z(izmin:izmax);
        timesTH(1).t=GridDan(1).Y1(:)'/1000;
        
        pdat(1).p=fact*squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[1]),3));
%                pdat(1).p=fact*squeeze(sum(TwoDDan(1).Q(izmin:izmax,:,[1]),3));

        imaxovr=[1 0];
        maxcovOvr=8;
        iminovr=[1 0];
        mincovOvr=1;
        
        case 2
        [izmin izmax]=findheight(Grid.Z,14e3,22e3);
        zz(2).z=GridDan(1).Z(izmin:izmax);
        timesTH(2).t=GridDan(1).Y1(:)'/1000;
        
        pdat(2).p=squeeze(TwoDDan(1).TH2(izmin:izmax,:));
        %maxcovOvr=100;
        end 