function []=prectimser(rad,time,year,month,day,hrs,mins,hrstartles,ia,ib,a,b,iecho,n,draw,isave)

clear xlabel ylabel

ptype(1).p='Max';
ptype(2).p='Mode';
ptype(3).p='Variance';
ptype(4).p='Median';
ptype(5).p='Mean';
%iecho=2;
%isave=1;

dtype(1).d='Echo';
dtype(2).d='PPI';
dtype(3).d='Capppi';
dtype(4).d='Pmax';
dtype(5).d='MaxEcho';

for it=1:length(draw)
    ityp=draw(it);
    if ityp==3
        add='^2';
    else
        add='';
    end
exname=strcat('c:/matlabr12/work/bauru/casestudy/bauruGraphs/',ptype(ityp).p,dtype(iecho).d,'BIT');

%time=[-300:300:15*3600];
%n=length(Surf2)+1;
%starttimeind=3; %place where data in np begins
%hrstart=9;

switch iecho
case 1
    ylab=strcat(ptype(ityp).p,' of 10dBZ radar echo tops (km',add,')');
case 2
    ylab=strcat(ptype(ityp).p,' of Precipitation Rate (mm/hr)',add);
case 3
    ylab=strcat(ptype(ityp).p,' of Radar Reflectivity at 3.5km (dBZ',add,')');
case 4
    n=n-1;
    ylab=strcat(ptype(ityp).p,' of Precipitation Rate (mm/hr)',add);
case 5
    ylab=strcat(ptype(ityp).p,' of Max Radar Reflectivity Echoes (dBZ',add,')');
end

xlab='Local Time';

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
gcf=figure('name',strcat('timseries of:',ylab),'position',posit);

for i=1:n
   
    %ydat(i).y=rad(i).means;
    %ydat(i).y=squeeze(mean(Surf(i).instant,1));
    switch ityp
    case 1
        ydat(i).y=rad(i).maxs;
    case 2
        ydat(i).y=rad(i).modes;
    case 3
        ydat(i).y=rad(i).vars;
    case 4
        ydat(i).y=rad(i).means;
    case 5
        ydat(i).y=rad(i).medians;
    end
    
    xdat(i).x=time(ia:ib);  %:length(ydat(i).y));
    xdat(i).x=datenum(2004,2,24,hrstartles+floor(xdat(i).x/3600),60*(xdat(i).x/3600-floor(xdat(i).x/3600)),0); 
end


labs(1).l='CCN = 720cm^{-3}';
labs(2).l='CCN = 240cm^{-3}';

labs(3).l='c1000-1200s';   %114
labs(4).l='c125km';    %121
labs(5).l='c250km';    %121
labs(6).l='c500km';    %121
labs(7).l='c2000-500m';

if iecho~=4
	labs(n).l='Radar data';
	
	%a=find(hrs>=hrstart);
	%a=a(1);
	
	xdat(n).x=datenum(year(a:b),month(a:b),day(a:b),hrs(a:b),mins(a:b),0);
	%xdat(n).x=xdat(n).x(3:end);
	%ydat(n).y=ydat(n).y;
	%ydat(n).y=squeeze(meany2(27:end));
	%ydat(n).y=rms(27:end);
	
	
end

for i=1:n
    maxy(i)=max(ydat(i).y);
end
maxy=max(maxy);

if iecho==3 | iecho==5
    ylim([0 maxy]);
end


plotXY(xdat,ydat,labs,40,2);
tims=[9:2:23];
ti=datenum(2004,2,24,tims,0,0);
set(gca,'xtick',[ti]);
datetick('x',15,'keepticks');

xlabel(xlab);
ylabel(ylab);

if isave==1
    set(gcf,'paperpositionmode','auto');
	print(gcf,'-dbitmap','-r350',exname);
    %print(gcf,'-dmeta',exname);
	close(gcf);
end

end