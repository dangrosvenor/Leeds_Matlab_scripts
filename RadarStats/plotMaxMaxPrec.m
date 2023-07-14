%plots time height plots
fsize=18;
isave=1;

clear xdat ydat

for  i=1:2
exname=strcat('c:/matlabr12/work/bauru/casestudy/bauruGraphs/liqProf',num2str(i));
%xdat(i).x=9+time./3600;
%xdat(i).x=datenum(2004,2,24,hrstartles+floor(xdat(i).x/3600),60*(xdat(i).x/3600-floor(xdat(i).x/3600)),0);
ydat(i).y=GridDan(1).Z/1000;
xdat(i).x=squeeze(max(prof(i).liq,[],2));
scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

gcf=figure('position',posit);
labs(2).l='CCN = 720cm^{-3}';
labs(1).l='CCN = 240cm^{-3}'; %profs were loaded in this order (240 then 720)
end

plotXY(xdat,ydat,labs,40);




% 
% tims=[9:2:23];
% ti=datenum(2004,2,24,tims,0,0);
% set(gca,'xtick',[ti]);
% datetick('x',15,'keepticks');


set(gca,'fontsize',fsize);
xlabel('Mean Liquid Water (kg/kg)');
ylabel('Height (km)');
%title('Percentage of Original Maximum Tracer (%)');

if isave==1
     set(gcf,'paperpositionmode','auto');
	print(gcf,'-dbitmap','-r350',exname);
    %print(gcf,'-dmeta',exname);
	%close(gcf);
end
