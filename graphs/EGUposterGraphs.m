clear labs xdat ydat

figname='EGU graph';
gridon=1; %switch for grid default =1 - note probs with grid when resizing due to extra ticks.s
scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.13];
zmax=25;
lor=0;
noplot=0;
ylogflag=0;


 Ms=1/6*pi*rhoS*Ds.^3;


a=find(times~=0);
dt(2:length(times))=times(2:end)-times(1:end-1);
dt(a(1))=dt(a(2));
a=find(times==0);
dt(a)=1; %avoid divide by zero
dtt=repmat(dt,[150 1]);
dtt2=repmat(dtt,[1 1 14]);
dtt3=permute(dtt2,[1 3 2]);





graph=2;

switch graph
    
    case 10
          
      
    sumUnDep=sum(NaeroDetQ(:,7:17,:)*500e3^2,2); %now (149,1,39)
    sumDep=sum(depN(:,7:17,:)*500e3^2,2); %now (149,1,39)]
    
    for iq=1:39
        
        a=find(sumUnDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosUnDep(iq)=sum(sumUnDep(a,1,iq),1);
        
        a=find(sumDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosDep(iq)=sum(sumDep(a,1,iq),1);
        
    end
    
    figname='Aerosol Number Proportion';
        %ydat(1).y=(sumPosDep-sumPosUnDep)./sumPosDep;    
    ydat(1).y=(sumPosUnDep)./sumPosDep;
    ydat(1).y=max([ydat(1).y; zeros([1 39])]);
    
    xdat(1).x=(Ms(2:end)+Ms(1:end-1)/2);
    xdat(1).x=Dsmid;
    
	labs(1).l='Ratio of detrained aerosol number in activated case to unactivated case';
	
	xlab='Aerosol diameter (m)';
    %xlab='Aerosol mass (kg)';
	ylab='Fraction';
    
    
    zmax=1.2;
    
    ylogflag=1;
    lor=2;
    gridon=0;
    
    %set(gca,'xlim',[9e-9 2e-6]); do this after


    
    case 9
        
    noplot=1;    
      
    sumUnDep=sum(NaeroDetQ(:,7:17,:)*500e3^2,2); %now (149,1,39)
    sumDep=sum(depN(:,7:17,:)*500e3^2,2); %now (149,1,39)]
    
    for iq=1:39
        
        a=find(sumUnDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosUnDep(iq)=sum(sumUnDep(a,1,iq),1);
        
        a=find(sumDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosDep(iq)=sum(sumDep(a,1,iq),1);
        
    end
    
    figname='Aerosol Number Histogram';
    
    figure('name',figname,'Position',posit);
   
     Ms=Ds;

	dm=-log(Ms(2:end))+log(Ms(1:end-1));
    
   
    


    dat=sumPosDep./dm;
    h=stairs((Ms(1:end-1)),dat);
    set(h,'linewidth',4,'color','b');
    
    hold on
    
    dat2=sumPosUnDep./dm;
    h=stairs((Ms(1:end-1)),dat2);
    set(h,'linewidth',4,'color','r');
    
    
   
     for i=1:39
        h2=line([(Ms(i)) (Ms(i))],[0 dat(i)]);
        set(h2,'linewidth',4);
        set(h2,'color','b');  
    end

    for i=1:39
        h=line([(Ms(i)) (Ms(i))],[0 dat2(i)]);
        set(h,'linewidth',4);
        set(h,'color','r');  
    end
   
    for i=2:40
        h=line([(Ms(i)) (Ms(i))],[0 dat2(i-1)]);
        set(h,'linewidth',4);
        set(h,'color','r');  
    end
   
	
	set(gca,'xscale','log');
    %set(gca,'ticklength',[0.005 0.005])
    set(gca,'fontsize',23)
    
	labs(1).l='With activation';
	labs(2).l='Without activation';
	
	xlab='Aerosol mass (kg)';
    xlab='Aerosol diameter (m)';
	ylab='Total detrained aerosol number divided by log bin width (#/m)';
    
    xlabel(xlab);
    ylabel(ylab);
    
    legend([h h2],labs.l,1);

    
    set(gca,'xlim',[9e-9 2e-6]);
    
    Ms=1/6*pi*rhoS*Ds.^3;
    
    
    case 8
          
      
    sumUnDep=sum(MaeroDetQ(:,7:17,:)*500e3^2,2); %now (149,1,39)
    sumDep=sum(depM(:,7:17,:)*500e3^2,2); %now (149,1,39)]
    
    for iq=1:39
        
        a=find(sumUnDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosUnDep(iq)=sum(sumUnDep(a,1,iq),1);
        
        a=find(sumDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosDep(iq)=sum(sumDep(a,1,iq),1);
        
    end
    
    figname='Aerosol Mass Proportion';
        %ydat(1).y=(sumPosDep-sumPosUnDep)./sumPosDep;    
    ydat(1).y=(sumPosUnDep)./sumPosDep;
    ydat(1).y=max([ydat(1).y; zeros([1 39])]);
    
    xdat(1).x=(Ms(2:end)+Ms(1:end-1)/2);
    xdat(1).x=Dsmid;
    
	labs(1).l='Ratio of detrained aerosol mass in activated case to unactivated case';
	
	xlab='Aerosol mass (kg)';
    xlab='Aerosol diameter (m)';
	ylab='Fraction';
    
    
    zmax=1.2;
    
    ylogflag=1;
    lor=2;
    gridon=0;

    %set(gca,'xlim',[9e-9 2e-6]) - do this after
    
    case 7
        
    noplot=1;    
      
    sumUnDep=sum(MaeroDetQ(:,7:17,:)*500e3^2,2); %now (149,1,39)
    sumDep=sum(depM(:,7:17,:)*500e3^2,2); %now (149,1,39)]
    
    for iq=1:39
        
        a=find(sumUnDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosUnDep(iq)=sum(sumUnDep(a,1,iq),1);
        
        a=find(sumDep(:,1,iq)>0);
        b=find(a>35);
        a=a(b); %only points above cloud base
        sumPosDep(iq)=sum(sumDep(a,1,iq),1);
        
    end
    
    figname='Aerosol Mass Histogram';
    
    figure('name',figname,'Position',posit);

    
        Ms=Ds;
	dm=-log(Ms(2:end))+log(Ms(1:end-1));
    
    


    dat=sumPosDep./dm;
    h=stairs((Ms(1:end-1)),dat);
    set(h,'linewidth',4,'color','b');
    
    hold on
    
    dat2=sumPosUnDep./dm;
    h=stairs((Ms(1:end-1)),dat2);
    set(h,'linewidth',4,'color','r');
   
     for i=1:39
        h2=line([(Ms(i)) (Ms(i))],[0 dat(i)]);
        set(h2,'linewidth',4);
        set(h2,'color','b');  
    end

    for i=1:39
        h=line([(Ms(i)) (Ms(i))],[0 dat2(i)]);
        set(h,'linewidth',4);
        set(h,'color','r');  
    end
   
    for i=2:40
        h=line([(Ms(i)) (Ms(i))],[0 dat2(i-1)]);
        set(h,'linewidth',4);
        set(h,'color','r');  
    end
   
	
	set(gca,'xscale','log');
    %set(gca,'ticklength',[0.005 0.005])
    set(gca,'fontsize',23)
    
	labs(1).l='With activation';
	labs(2).l='Without activation';
	
	xlab='Aerosol mass (kg)';
    xlab='Aerosol diameter (m)';
	ylab='Total detrained aerosol mass divided by log bin width (kg/m)';
    
    xlabel(xlab);
    ylabel(ylab);
    
    legend([h h2],labs.l,2);
    
    zmax=17.5;
    
    set(gca,'xlim',[9e-9 2e-6]);
    
    Ms=Msbk;
    
    
    case 6
        
    dz=GridDan(1).Z(2:end)-GridDan(1).Z(1:end-1);
    dzz=repmat(dz,[1 size(diag(1).dg,3) size(TwoDDan(1).Q,3)-15]);
        
    sumMaero=sum(MaeroDetQ*500e3^2./dzz(:,:,:),3);
    cumMaero=cumsum(sumMaero(:,7:17),2);
    
    sumUnDep=sum(depM*500e3^2./dzz(:,:,:),3);
    cumUnDep=cumsum(sumUnDep(:,7:17),2);
    
	xdat(1).x=cumMaero(:,11);
	ydat(1).y=Grid.Z(2:end)./1000;
	
	xdat(2).x=cumUnDep(:,11);
	ydat(2).y=Grid.Z(2:end)./1000;
	
	
	labs(1).l='With activation';
	labs(2).l='Without activation';
	
	xlab='Loss of aerosol mass from vertical mass flux (kg/m)';
	ylab='Height (km)';
    
    zmax=17.5;
    
    
    case 5
    dz=GridDan(1).Z(2:end)-GridDan(1).Z(1:end-1);
    dzz=repmat(dz,[1 size(diag(1).dg,3) size(TwoDDan(1).Q,3)]);
    dzz=permute(dzz,[1 3 2]);

    summ=cumsum(MDetQ(1).m(:,:,50:68),3)*500e3^2./dzz(:,:,50:68);
    summ2=cumsum(MDetQ(2).m(:,:,50:68),3)*500e3^2./dzz(:,:,50:68);
    
    summ=sum(summ(:,2:6,:),2);
    summ2=sum(summ2(:,2:6,:),2);
    
     
	xdat(1).x=summ(:,end);
	ydat(1).y=Grid.Z(2:end)./1000;
	
	xdat(2).x=summ2(:,end);
	ydat(2).y=Grid.Z(2:end)./1000;
	
	
	labs(1).l='Low updraught case';
	labs(2).l='High updraught case';
	
    
	xlab='Detrainment of condensate mass per unit height (kg/m)';
    
	ylab='Height (km)';
    
    figname='Condensate Mass Gain';
    
    zmax=19;
    lor=0;
   
    
    case 4
        
        dz=GridDan(1).Z(2:end)-GridDan(1).Z(1:end-1);
		dzz=repmat(dz,[1 size(diag(1).dg,3) size(TwoDDan(1).Q,3)]);
        dzz=permute(dzz,[1 3 2]);
    
    summ=cumsum(MDetQ(1).m(:,:,50:68),3)*500e3^2./dzz(:,:,50:68);
    summ2=cumsum(MDetQ(2).m(:,:,50:68),3)*500e3^2./dzz(:,:,50:68);
    
     
	xdat(1).x=summ(:,1,end);
	ydat(1).y=Grid.Z(2:end)./1000;
	
	xdat(2).x=summ2(:,1,end);
	ydat(2).y=Grid.Z(2:end)./1000;
	
	
	labs(1).l='Low updraught case';
	labs(2).l='High updraught case';
	
    
	xlab='Detrainment of water vapour mass per unit height (kg/m)';
    
	ylab='Height (km)';
    
    figname='Vapour Mass Gain';
    
    zmax=19;
    lor=0;
    
    
case 3
    
    summ=cumsum(flux(1).m(:,:,50:68)./dtt3(:,:,50:68),3);
    summ2=cumsum(flux(2).m(:,:,50:68)./dtt3(:,:,50:68),3);
    
    
	xdat(1).x=summ(:,14,end);
	ydat(1).y=Grid.Z(1:end)./1000;
	
	xdat(2).x=summ2(:,14,end);
	ydat(2).y=Grid.Z(1:end)./1000;
	
	
	labs(1).l='Low updraught case';
	labs(2).l='High updraught case';
	
    
    xlab='Average upwards ozone flux (kg/m^2/s)';
    
	ylab='Height (km)';
    
    figname='Ozone Flux';
    
case 2
    
    dz=GridDan(1).Z(2:end)-GridDan(1).Z(1:end-1);
	dzz=repmat(dz,[1 size(diag(1).dg,3) size(TwoDDan(1).Q,3)]);
    dzz=permute(dzz,[1 3 2]);
    
    summ=cumsum(MDetQ(1).m(:,:,50:68),3)*500e3^2./dzz(:,:,50:68);
    summ2=cumsum(MDetQ(2).m(:,:,50:68),3)*500e3^2./dzz(:,:,50:68);
    
     
	xdat(1).x=summ(:,14,end);
	ydat(1).y=Grid.Z(2:end)./1000;
	
	xdat(2).x=summ2(:,14,end);
	ydat(2).y=Grid.Z(2:end)./1000;
	
	
	labs(1).l='Low updraught case';
	labs(2).l='High updraught case';
	
    
	xlab='Detrainment of ozone mass per unit height (kg/m)';
    
	ylab='Height (km)';
    
    figname='Ozone Mass Gain';
    
    zmax=19;
    lor=4;
   

case 1
	xdat(1).x=summ(2:end)*300.*dz;
	ydat(1).y=Grid.Z(2:end)./1000;
	
	xdat(2).x=sumMdet;
	ydat(2).y=Grid.Z(2:end)./1000;
	
	
	labs(1).l='Activation loss';
	labs(2).l='Increase w/ height due to mass flux';
	
	xlab='Increase in Aerosol Number Flux w/ Height (#/m^2)';
	xlab='Aerosol Mass (kg/m^2)';
	ylab='Height (km)';
end


if noplot==0
figure('name',figname,'Position',posit);
[h,ax,ax2]=plotXY3(xdat,ydat,labs,0,4,lor,ylogflag,xlab,ylab,[0 zmax],1);
%function [H1,ax2]=plotXY3(xdat,ydat,labs,nmark,lwidth,leglor,logflag,xlab,ylab,ylims,zline) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers

% LEGEND(...,Pos) places the legend in the specified
%     location:
%         0 = Automatic "best" placement (least conflict with data)
%         1 = Upper right-hand corner (default)
%         2 = Upper left-hand corner
%         3 = Lower left-hand corner
%         4 = Lower right-hand corner
%        -1 = To the right of the plot

if gridon==1
    grid on;
end

end

