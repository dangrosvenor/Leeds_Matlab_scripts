%plots several different lines on one graph 17/10/03 14:10
figure('name','Combination plot');
titsbaurufull;


conm(1)=61;
conm(2)=62;
conm(3)=23;

size(conm);
sconm=ans(2);
pnames;
%nplots=input('Enter number of plots :');
%con=input('Enter col no.: ');
%fact=input('Please enter scale factor: ');
nplots=12;

scrsz=get(0,'ScreenSize');

clear minsx maxsx minsy maxsy;
for i=1:nplots
    minsx(i)=min(SerDan(i).SER(:,1));
    maxsx(i)=max(SerDan(i).SER(:,1));
    minsy(i)=min(SerDan(i).SER(:,con));
    maxsy(i)=max(SerDan(i).SER(:,con));    
end;
minx=min(minsx);
maxx=max(maxsx);
miny=min(minsy);
maxy=max(maxsy)*1.1;


exist posit;
a=ans;
if a==0;
    posit=[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2.2];
end;

for k=1:1
    cons=int2str(con);
    %figure('name',strcat('Column:',cons,' v.time'),'Position',posit);
    
    jmax=6;
    if nplots-k*6<0 
        jmax=rem(nplots,6)+(k-1)*6;
    end;
    
    for j=1:12;
        plotnum=j;
        if plotnum==0;
            plotnum=6;
        end
    hold on;
    subplot(6,2,plotnum);
    hold on;
    
    hs(j).h=newplot;
    
    times=SerDan(j).SER(:,1);
    mlt=SerDan(j).SER(:,61)+SerDan(j).SER(:,66)+SerDan(j).SER(:,65)+SerDan(j).SER(:,63); %graupel, ice and snow melting + graupel shedding
    
    [AX,H1,H2] = plotyy(times,mlt,times,SerDan(j).SER(:,conm(3)));
    set(H1,'LineStyle','--');
    %set(H2,'LineStyle',':');
    
    H3=plot(SerDan(j).SER(:,1),SerDan(j).SER(:,conm(2)),'r:');
    H4=plot(SerDan(j).SER(:,1),SerDan(j).SER(:,64),'k-.');
    if j==11
        set(get(AX(1),'Ylabel'),'String','Total Process Rate (kg/m^2/s)');
        set(get(AX(2),'Ylabel'),'String','Total Rain Mixing Ratio (kg/m^2)');
        lh=legend([H1 H3 H2],'Melting and graupel shedding process rates','Autoconversion process rate','Total rain mixing ratio',2);
        set(lh,'Fontsize',7);
    end
    
    set(AX,'Fontsize',7);
    set(AX,'xlim',[0 maxx]);
    %axis([minx maxx miny maxy]);
    title(tit(j).tit,'Fontsize',8);
    xlabel('Time (s)','Fontsize',7,'position',[5000 -0.8e-4]);
    
    %ylabel('Precipitation (mm/hour)');
    %ylabel(strcat('Column ',cons));
    end
end

%text(-8,21.2,textdataDan(1).text,'units','centimeters');

for j=1:nplots
    dat=SerDan(j).SER;
    TrapIntDatN2;
end;
if con>=60
    figlab=strcat('Averages of Time Series for: ',pname(con).p,',',textdataDan(1).text);
else
    figlab=strcat('Averages of Time Series for column: ',cons,',',textdataDan(1).text);
end

%ccnplN;
