updown=1;

    logflag=0;
    
    
    
    figname='Sources';
    
    xlims=1;
    xlimits=[-20 20];
    
    zmin=12;  %(km)
    zmax=30;   %z(end)/1000;
    

    z=GridDan(idir).Z;
    
    
    dz=GridDan(idir).Z(2:end)-GridDan(idir).Z(1:end-1); 
    rho=GridDan(idir).RHO;
    
    t1=3;
    t2=86; %index for start and end times of difference plot
        iz=1;
        iz2=length(z);
        
        switch updown
        case 1
        
        xdat(1).x=f*( sum(icediag4(idir).i(:,t2,[35]),3) - sum(icediag4(idir).i(:,t1,[35]),3) ); %difference in end and start vapour diags
        xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediag4(idir).i(:,t1:t2,27),3),GridDan(idir).t,iz,iz2); %microphys
        xdat(3).x=xdat(1).x-xdat(2).x;
        xdat(4).x=f*TotMassBudgetProf(GridDan(idir),sum(icediag4(idir).i(:,t1:t2,[1]),3),GridDan(idir).t,iz+1,iz2); %microphys
        
        xlab='Source of Vapour from Updraughts (ppmv)';
        
        case 2
        xdat(1).x=f*( sum(icediag4(idir).i(:,t2,[36]),3) - sum(icediag4(idir).i(:,t1,[36]),3) ); %difference in end and start vapour diags
        xdat(2).x=f*TotMassBudgetProfRate2(GridDan(idir),sum(icediag4(idir).i(:,:,28),3),GridDan(idir).t,iz,iz2);
        xdat(3).x=xdat(1).x-xdat(2).x;
        
        xlab='Source of Vapour from Downdraughts (ppmv)';
        
        end
        
		ylab='Height (km)';

        
       
    
    ydat(1).y=z/1000;
    labs(1).l='Vapour Change';
    
    ydat(2).y=z/1000;
    labs(2).l='Microphysics';
    
    ydat(3).y=z/1000;
    labs(3).l='Advection';
    
    ydat(4).y=z(2:end-1)/1000;
    labs(4).l='Advection';
    
%     ydat(3).y=Grid.Z(2:len+1)/1000;
%     labs(3).l='Microphysics';
%     
%     ydat(4).y=Grid.Z(2:len+1)/1000;
%     labs(4).l='Fall Speed Flux';
%     
%     ydat(5).y=Grid.Z(2:len+1)/1000;
%     labs(5).l='Ice Flux';
%     
%     ydat(6).y=Grid.Z(2:len+1)/1000;
%     labs(6).l='Upwards Ice Flux';
%     
%     ydat(7).y=Grid.Z(2:len+1)/1000;
%     labs(7).l='Downwards Ice Flux';
%     
%     ydat(8).y=Grid.Z(2:len+1)/1000;
%     labs(8).l='Net Vapour Flux';
%     
%     ydat(9).y=Grid.Z(2:len+1)/1000;
%     labs(9).l='Ice Fall + Flux';
    
    %xdat(2).x =f*( -m(1).m + m(11).m(2:end) ); %vapour gained from detrainment + microphysical sources
                                   %calculated with VapBudget.m
                                   
%     [iz,iz2]=findheight(Grid.Z,0e3,30e3);
%     ydat(2).y=Grid.Z(iz+1:iz2)/1000;
%     labs(2).l='Net Vapour Gain Flux - Microphysics';
    
    
    lor=1;
    
    %stuff for additional pressure axis
    secyA=z/1000;
    secyB=GridDan(idir).PREFN/100;
    lab2='Pressure (hPa)';  
    dual=1;