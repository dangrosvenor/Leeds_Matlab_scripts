%idir=1;
    
	logflag=0;
    fact=1e6*28.97/18;
    
    %iminovr=1;
%    mincovOvr=(1e-7);
    mincovOvr=(-7);

    
    %imaxovr=1;
    maxcovOvr=(0.05);
    maxcovOvr=(5E-2);

        
	tit(1).tit='Mean Fall Speed Flux (m/s)';
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
    
  
  
   
 
   hm='fall ratio mr';
   
   switch hm
       
   case 'fall ratio mr'
       iceloss=f*TotMassBudgetALL(GridDan(idir),sum(icediag4(idir).i(:,3:86,[11 12 17 18 23 24]),3),GridDan(idir).t,izmin,izmax+1);
        
       %pdat(i).p=f*sum(icediag4(idir).i(izmin:izmax,3:86,[37:42]),3) - iceloss;
       
       pdat(i).p=iceloss*3600/300; %ice loss in ppmv/hr
       
       i0=find(pdat(i).p>0);
       pdat(i).p(i0)=1e-10;
       pdat(i).p=abs(pdat(i).p);
       %ipos=find(pdat(1).p>0);
       %pdat(i).p(ipos)=0;
       
      % tit(1).tit='Total Water After Fall Speed Removal';
       tit(1).tit='Fall Speed Removal Rate (ppmv/hr)';

       figlab='Fall speed removal rate';
       
        iminovr=1;
        mincovOvr=(-3);
        
        %imaxovr=0;
        maxcovOvr=(0.08);
        
        ncont=20;
        
        sig=2;
        
        clines=0;
        clab=1;
        
        logflag=1;
       
   case 'fall ratio'
%        pdat(i).p=TotMassBudgetALL(GridDan(idir),sum(icediag4(idir).i(:,3:86,[11 12 17 18 23 24]),3),GridDan(idir).t,izmin,izmax+1) ...
%            ./ sum(icediag4(idir).i(izmin:izmax,3:86,[37:42]),3);

         pdat(i).p=fact*3600/300*TotMassBudgetALL(GridDan(idir),sum(icediag4(idir).i(:,3:86,[11 12 17 18 23 24]),3),GridDan(idir).t,izmin,izmax+1);
        
       %ipos=find(pdat(1).p>0);
       %pdat(i).p(ipos)=0;
       
       tit(1).tit='Fall Speed Removal Fraction';
       figlab='Fall Ratio time-height';
       
       tit(1).tit='Fall Speed Removal Rate (ppmv/hr)';
       figlab='Fall speed removal time-height';
       
       
        iminovr=1;
        mincovOvr=(-200);
        
        %imaxovr=1;
        maxcovOvr=(0);
        
        ncont=20;
        
        sig=1;
        
        clines=0;
        clab=0;
        
        sig=2;
       
   case 'mean flux'
       pdat(i).p=fact*squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,[11 12 17 18 23 24]),3));
       
   case 'all'
      
    pdat(i).p=squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,[11 12 17 18 23 24]),3) ...
        ./ sum(icediag4(idir).i(izmin:izmax,dumprange,[37:42]),3) );
    
	case 'snow'

    pdat(i).p=squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,[11 12]),3) ...
        ./ sum(icediag4(idir).i(izmin:izmax,dumprange,[37:38]),3) );
    
	case 'graupel'
    
    pdat(i).p=squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,[17 18]),3) ...
        ./ sum(icediag4(idir).i(izmin:izmax,dumprange,[39:40]),3) );
    
	case 'ice'
    pdat(i).p=squeeze(sum(icediag4(idir).i(izmin:izmax,dumprange,[23 24]),3) ...
        ./ sum(icediag4(idir).i(izmin:izmax,dumprange,[41:42]),3) );
    
end