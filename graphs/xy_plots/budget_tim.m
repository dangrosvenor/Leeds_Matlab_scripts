
    logflag=0;
    
    dumprange=[2:87];
    %z=GridDan(idir).Z;
 
    xlab=['Time (UTC)'];
	ylab='Total Change from Top Down to 16km';
    
    

    xdat(1).x = GridDan(idir).t(dumprange)+3;
    labs(1).l='Microphysical Vapour Source';
    
    xdat(2).x= GridDan(idir).t(dumprange)+3;
    labs(2).l='Fall Speed Flux';
    
    xdat(3).x= GridDan(idir).t(dumprange)+3;
    labs(3).l='Overall Vapour Change';
    
    xdat(4).x= GridDan(idir).t(dumprange)+3;
    labs(4).l='Overall Tot Water Change';
    
%     
%     secyA=z/1000 + 0.62;
%     secyB=GridDan(idir).PREFN/100;
%     lab2='Pressure (hPa)';  
%     dual=1;
%    

    zbttl=16e3;
    figname=['top down budget to ' num2str(zbttl/1000) 'km ' direc(idir).dir];

    %[za zb]=findheight(GridDan(idir).Z,16e3-620,16.5e3-620);
    [za]=findheight(GridDan(idir).Z,zbttl-620);
    zb=length(GridDan(idir).Z);
    %zb=za;
    
     izlim=0; %flag to say that don't want to scale z axis   
     zmin=0;  %(km)
     zmax=0;
     
     rho=GridDan(idir).RHON;
     dz=GridDan(idir).Z(2:end)-GridDan(idir).Z(1:end-1);
     
     icum=0;
     
    if idir==2 
        for it=dumprange %run through dumprange
            ydat(1).y(it-dumprange(1)+1)=sum( icediagsALL(idir).i(za:zb,it,28).*rho(za-1:zb-1).*dz(za-1:zb-1) ,1 ) ;
            ydat(2).y(it-dumprange(1)+1)=sum( TotMassBudgetProf(GridDan(idir),sum(icediagsALL(idir).i(:,it,22:24),3),GridDan(idir).t,za,zb) ,1 ) ;
            ydat(3).y(it-dumprange(1)+1)=sum( (icediagsALL(idir).i(za:zb,it+1,37) - icediagsALL(idir).i(za:zb,it-1,37) )...
                .*rho(za-1:zb-1).*dz(za-1:zb-1) , 1 )/600 ; %vapour
            
            ydat(4).y(it-dumprange(1)+1)=sum( ( sum(icediagsALL(idir).i(za:zb,it+1,[10:15 37:42]),3) - sum(icediagsALL(idir).i(za:zb,it-1,[10:15 37:42]),3) )...
                .*rho(za-1:zb-1).*dz(za-1:zb-1) , 1 )/600 ; %vapour
            
        end
        
        iylen=icum*length(ydat);
        for iy=1:iylen
            ydat(iy).y=cumsum(ydat(iy).y);
        end
        
    end
%     
    
     xlims=1;
     xlimits=[xdat(1).x(1) xdat(1).x(end)];     
     %nmark=-1;