%plot distributions

isave=1;
exname='c:/matlabr12/work/bauru/casestudy/bauruGraphs/CappiDist_14-23BIT';

n=length(Surf2)+1;

if iecho==2
    ylab='Log10 of Normalised Distribution';
    xlab='Precipitation Rate (mm/hr)';
elseif iecho==1
    ylab='Log10 of Normalised Distribution of Echo Tops';
    xlab='10dBZ Radar Echo Top (km)';
elseif iecho==3
    ylab='Log10 of Normalised Distribution';
    xlab='Radar Reflectivity at 3.5km (dBZ)';
else
    ylab='Log10 of Normalised Distribution';
    xlab='Radar Reflectivity of Vertical Radar Slices (dBZ)';
end



%time=[-300:300:15*3600];
%imported 1-169 so = 3-170; pos 3 = 300secs

onetime=0; %flag for plotting distribution at only one time 
%hrstartles=9; %start time for les data
timday=17; %time of day for plot

tim=timday-hrstartles; %no. hours after hrstart=time(1) for plot




% a=find(hrs>=hrstart); %nps data starts at hrstart time so need to subtract index of this
% a=a(1);

b=find(hrs>=timday); %for radar data
b=b(1)-a+1;

[temp itim]=min( abs((tim)*3600 - time) ); %for LES data
itim=itim-ia+1;

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
gcf=figure('name','Prec Size Distribution','position',posit);
for i=1:n
    xdat(i).x=rad(i).mid;
    
    if onetime==0
        ydat(i).y=log10(rad(i).np./rad(i).tot./rad(i).diffs);
    else
        if i~=n
            ydat(i).y=log10(rad(i).nps(itim,:)./sum(rad(i).nps(itim,:))./rad(i).diffs);
        else
            ydat(i).y=log10(rad(i).nps(b,:)./sum(rad(i).nps(b,:))./rad(i).diffs);
        end
    end
    
     if iecho==1
         xdat(i).x=xdat(i).x(1:end-3); %because have 2 values at end relating to unknown heights
         ydat(i).y=ydat(i).y(1:end-3);
     end
end



labs(1).l='CCN = 720cm^{-3}';
labs(2).l='CCN = 240cm^{-3}';
labs(3).l='c1000-1200s';   %114
labs(4).l='c125km';    %121
labs(5).l='c250km';    %121
labs(6).l='c500km';    %121
labs(7).l='c2000-500m'; 
labs(n).l='Radar data';


plotXY(xdat,ydat,labs,40);

xlabel(xlab);
ylabel(ylab);

if isave==1
    set(gcf,'paperpositionmode','auto');
	%print(gcf,'-djpeg','-r350',exname);
    print(gcf,'-dbitmap',exname);
	%close(gcf);
end