scrsz=get(0,'ScreenSize');
posit=[9 214 scrsz(3)/1.01 scrsz(4)/1.4];
%figure('name','actGraphs','Position',posit);

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
figure('name','actGraphs','Position',posit);

graph=2;


switch graph
    
case 2
    sumMaero=sum(MaeroDetQ,3);
    cumMaero=cumsum(sumMaero(:,7:17),2);
    
    sumUnDep=sum(depM,3);
    cumUnDep=cumsum(sumUnDep(:,7:17),2);
    
	xdat(1).x=cumMaero(:,11);
	ydat(1).y=Grid.Z(2:end)./1000;
	
	xdat(2).x=cumUnDep(:,11);
	ydat(2).y=Grid.Z(2:end)./1000;
	
	
	labs(1).l='With activation';
	labs(2).l='Without activation';
	
	xlab='Increase in aerosol number flux w/ height (#/m^2)';
	xlab='Gain in aerosol mass flux w/ height (kg/m^2)';
	ylab='Height (km)';
   


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

plotXY3(xdat,ydat,labs,0,4,2,0,xlab,ylab,[0 19]);
grid on;