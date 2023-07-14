%plots time height plots
fsize=18;
isave=0;
%load c:/matlabr12/work/bauru/casestudy/forcecons/diag/profnew+ccn_2-169;

maxZ=18500;
maxtr=1.0;
time=[46*300:300:79*300];
hrstartles=9;


tit(1).tit='Max of Low Level Tracer';
tit(2).tit='Max of Low Level Tracer';

for  i=1:2
exname=strcat('c:/matlabr12/work/bauru/tracersjan2005/force+3_3th3qv/TracerTimH-',num2str(i));
%xdat(i).x=time;
%xdat(i).x=datenum(2004,2,24,hrstartles+floor(xdat(i).x/3600),60*(xdat(i).x/3600-floor(xdat(i).x/3600)),0);

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

iz=find(Grid.Z>=maxZ);
iz=iz(1);

gcf=figure('position',posit);
%pcolor(9+time./3600,Grid.Z(1:iz)./1000,maxLowTracer(i).prof(1:iz,47:80));hc=colorbar;%shading interp

pcolor(47:80,Grid.Z(1:iz)./1000,maxLowTracer(i).prof(1:iz,47:80));hc=colorbar;%shading interp
set(hc,'fontsize',fsize);




% 
% tims=[9:2:23];
% ti=datenum(2004,2,24,tims,0,0);
% set(gca,'xtick',[ti]);
% datetick('x',15,'keepticks');


set(gca,'fontsize',fsize);
xlabel('Local Time (hrs)');
ylabel('Height (km)');

title(tit(i).tit);

if isave==1
     set(gcf,'paperpositionmode','auto');
	print(gcf,'-djpeg','-r350',exname);
    %print(gcf,'-dmeta',exname);
	%close(gcf);
end

end