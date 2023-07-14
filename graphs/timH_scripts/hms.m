hm='fall';



%    minZ=12e3;
    minZ=14.7e3;
    maxZ=22e3;
    
    sig=2;
    clines=0; %makes black contour lines appear
    clab=1;
    
    ncont=25;
    
    %i2d=2; %tells it to label x axis in km
    
    z=GridDan(idir).Z; %change z for the different cases with kkp=230 for 25km and =250 for 30km tops
    time=GridDan(idir).t+3;
    
    
    
dumprange=[1:88];    
	logflag=0;
    fact=1e6*28.97/18;
    
    
    
    switch hm
    case 'vap'
	    tit(1).tit='Mean Vapour (ppmv)';
        figlab='Vapour time-height';
        pdat(i).p=f*sum(icediagALL(idir).i(izmin:izmax,dumprange,37),3);
        
        iminovr=1;
        mincovOvr=4;
        
        imaxovr=1;
        maxcovOvr=6;
         
    case 'total'
         tit(1).tit='Mean Total Water (ppmv)';
         figlab='Tot water time-height';
         pdat(i).p=f*sum(icediagALL(idir).i(izmin:izmax,dumprange,[37:42]),3);
         
         iminovr=1;
         mincovOvr=4;
        
        imaxovr=1;
        maxcovOvr=6;
         
     case 'ievap'
         tit(1).tit='Microphysical Source of Vapour (ppmv/hr)';
         figlab='Microphysical Source of Vapour';
%         pdat(i).p=cumsum(f*3600*sum(icediagALL(idir).i(izmin:izmax,dumprange,[28]),3),2);
          pdat(i).p=f*3600*sum(icediagALL(idir).i(izmin:izmax,dumprange,[28]),3);

         
         iminovr=0;
         mincovOvr=4;
        
        imaxovr=0;
        maxcovOvr=6;
        
     case 'fall'

       iceloss=f*TotMassBudgetALL(GridDan(idir),sum(icediagALL(idir).i(:,dumprange,[22 23 24]),3),GridDan(idir).t,izmin,izmax+1);
        
       %pdat(i).p=f*sum(icediag4(idir).i(izmin:izmax,3:86,[37:42]),3) - iceloss;
       
       pdat(i).p=iceloss*3600/300; %ice loss in ppmv/hr
       
       i0=find(pdat(i).p<0);
       pdat(i).p(i0)=1e-10;
       pdat(i).p=abs(pdat(i).p); %values of iceloss are negative for loss of ice - remove gain values for log plot
       %ipos=find(pdat(1).p>0);
       %pdat(i).p(ipos)=0;
       
      % tit(1).tit='Total Water After Fall Speed Removal';
       tit(1).tit='Fall Speed Removal Rate (ppmv/hr)';
       tit(1).tit='Fall Speed Source Rate (ppmv/hr)';
        
       figlab='Fall speed removal rate';
       
        iminovr=0;
        mincovOvr=(-3);
        
        imaxovr=0;
        maxcovOvr=(100);
        
        ncont=20;
        
        sig=2;
        
        clines=0;
        clab=0;
        
        logflag=0;
         
     end
  %  tit(1).tit='Min of Total Water (ppmv)';
    %tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	%tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    
    
