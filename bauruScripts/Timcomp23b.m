%plots graphs from N direcs (NO diff graphs & six to a page)

%nplots=input('Enter number of plots :');
%con=input('Enter col no.: ');

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

for k=1:fix((nplots-0.1)/6)+1
    cons=int2str(con);
    %figure('name',strcat('Column:',cons,' v.time'),'Position',posit);
    
    jmax=6;
    if nplots-k*6<0 
        jmax=rem(nplots,6)+(k-1)*6;
    end;
    
    for j=(k-1)*6+1:jmax;
        plotnum=rem(j,6);
        if plotnum==0;
            plotnum=6;
        end
    
    hs(j).h=subplot(rowtot,nplots,(rowno-1)*6+plotnum);
    %hs(j).h=newplot;
    size(direcDan(j).dir);
    ccns(j)=str2num(direcDan(j).dir(14:ans(2)));
    
    plot(SerDan(j).SER(:,1),SerDan(j).SER(:,con));
    axis([minx maxx miny maxy]);
    
    title(tit(j).tit);
    if rowno==3;
        xlabel('Time (s)');
    end;
    
    if j==1
    %ylabel(strcat('Column:',int2str(con)));    
    ylabel(ytit);
    end;
    set(gca,'Fontsize',7);
    %ylabel(strcat('Column ',cons));
    end
end

%text(-8,21.2,textdataDan(1).text,'units','centimeters');

for j=1:nplots
    dat=SerDan(j).SER;
    TrapIntDatN2;
end;
figlab=strcat('Averages of Time Series for column: ',cons);
%ccnplN;
