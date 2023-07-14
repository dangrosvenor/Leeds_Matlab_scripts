    logflag=0;
    
    z=GridDan(idir).Z/1000 +0.62;
 
    xlab=['Ice Sat MR (ppmv)'];
	ylab='Height (km)';
    
    figname=['ice sat mr - ' direc(idir).dir];
    

    ydat(1).y=z; %620m already added on in fourevents2.m
    labs(1).l='19.75 UTC';
 
    %T=tempLES(GridDan(idir)); %K
    %T=Tstart./(1e5./GridDan(idir).PREFN).^0.286;
    
    
%    dgfind=findhead('ALL_TH',dgstrDan(idir).dg);
% 	T=TimeAvDan(idir).DGAV(:,dgfind(1));
% 	T=T./(1e5./GridDan(idir).PREFN).^0.286;
%     xdat(1).x=SatVapPress(T,'goff','ice',GridDan(idir).PREFN,1); %Pa

    T=Tstart;
    xdat(1).x=SatVapPress(T,'goff','ice',GridDan(idir).PREFN,1); %Pa
    
    ydat(2).y=z; %620m already added on in fourevents2.m
    labs(2).l='2.8 UTC';
    
%     dgfind=findhead('ALL_TH',dgstrDan(idir).dg);
% 	T=TimeAvDan(idir).DGAV(:,dgfind(1));
% 	T=T./(1e5./GridDan(idir).PREFN).^0.286;
    
    T=Tend;
	xdat(2).x=SatVapPress(T,'goff','ice',GridDan(idir).PREFN,1); %Pa

    
    
    
% %     ydat(2).y=z(len:length(z)-1)/1000;
% %     labs(2).l='Vapour';
%     
     secyA=z;
     secyB=GridDan(idir).PREFN/100;
     lab2='Pressure (hPa)';  
     dual=1;
%     
     %xlims=1;
     xlimits=[0 10];
%     
     zmin=14;  %(km)
     zmax=19;
     
     %nmark=-1;



