function []=prectimser(rad,timestart,timeend,isave,draw)

clear xlabel ylabel

ptype(1).p='Max';
ptype(2).p='Mode';
ptype(3).p='Variance';
ptype(4).p='Median';
ptype(5).p='Mean';
ptype(6).p='Max Echo Tops';
%iecho=2;
%isave=1;

% dtype(1).d='Echo';
% dtype(2).d='PPI';
% dtype(3).d='Capppi';
% dtype(4).d='Pmax';
% dtype(5).d='MaxEcho';

for it=1:length(draw)  %draw contains the types of graphs to draw, e.g. =[1 2 3 4 5] for max,mode,etc.
    ityp=draw(it);
    if ityp==3
        add='^2';
    else
        add='';
    end
    
exname=strcat('c:/matlabr12/work/bauru/casestudy/bauruGraphs/slices/',ptype(ityp).p,'_echoTops');

%time=[-300:300:15*3600];
%n=length(Surf2)+1;
%starttimeind=3; %place where data in np begins
%hrstart=9;

ylab=strcat(ptype(ityp).p,' of Radar Slice Reflecivities (dBZ',add,')');
xlab='Local Time';

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
gcf=figure('name',strcat('timseries of:',ylab),'position',posit);

n=length(rad);
for i=1:n
   
    a=find(rad(i).hrs>=timestart & rad(i).hrs<timeend);
    b=a(end);
    a=a(1);
    
    %ydat(i).y=rad(i).means;
    %ydat(i).y=squeeze(mean(Surf(i).instant,1));
    switch ityp
    case 1
        ydat(i).y=rad(i).maxs(a:b);
    case 2
        ydat(i).y=rad(i).modes(a:b);
    case 3
        ydat(i).y=rad(i).vars(a:b);
    case 4
        ydat(i).y=rad(i).means(a:b);
    case 5
        ydat(i).y=rad(i).medians(a:b);
    case 6
        ylab=strcat(ptype(ityp).p,' of Radar Slice Reflecivities (km',add,')');
        ydat(i).y=rad(i).echotops(a:b);
    case 7
        ydat(i).y=rad(i).maxs2(a:b);
    end
    
    
    xdat(i).x=datenum(2004,2,24,rad(i).hrs(a:b),rad(i).mins(a:b),0);
end

    

labs(1).l='CCN = 240cm^{-3}';
labs(2).l='CCN = 720cm^{-3}';

labs(3).l='c1000-1200s';   %114
labs(4).l='c125km';    %121
labs(5).l='c250km';    %121
labs(6).l='c500km';    %121
labs(7).l='c2000-500m';

for i=3:7
	labs(i).l=strcat('Radar slice ',int2str(i-2));
end
	
	
	

for i=1:n
    maxy(i)=max(ydat(i).y);
end
maxy=max(maxy);
ylim([0 maxy]);
    


plotXY(xdat,ydat,labs,40,2);
tims=[timestart:2:timeend];
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

end %for it=1:length(draw)