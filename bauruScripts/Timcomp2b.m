%plots graphs from N direcs (NO diff graphs & six to a page)

%nplots=input('Enter number of plots :');
%con=input('Enter col no.: ');

scrsz=get(0,'ScreenSize');
clear dat;


clear minsx maxsx minsy maxsy;
for i=nstart:nend
    if con==99
        dat=SerDan(i).SER(:,18)+SerDan(i).SER(:,20)+SerDan(i).SER(:,21);
    else
        dat=SerDan(i).SER(:,con);
    end;
    minsx(i)=min(SerDan(i).SER(:,1)./3600);
    maxsx(i)=max(SerDan(i).SER(:,1)./3600);
    minsy(i)=min(dat);
    maxsy(i)=max(dat);    
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
    
    jmax=12;
    if nplots-k*6<0 
        jmax=rem(nplots,6)+(k-1)*6;
    end;
    
    for j=nstart:nend
      
        %plotnum=rem(j,6);
        plotnum=j;
        if plotnum==0;
            plotnum=6;
        end
    
    hs(j).h=subplot(rowtot,nplots,(rowno-1)*jmax+plotnum);
    %hs(j).h=newplot;
    size(direcDan(j).dir);
    ccns(j)=str2num(direcDan(j).dir(43:ans(2)));
    
    
    if con==99
        dat=SerDan(j).SER(:,18)+SerDan(j).SER(:,20)+SerDan(j).SER(:,21);
    else
        dat=SerDan(j).SER(:,con);
    end;
    
    plot(SerDan(j).SER(:,1)./3600,dat);
    axis([minx maxx miny maxy]);
    
    title(tit(j).tit,'Fontsize',5);
    if rowno==rowtot;
        xlabel('Time (hrs)','Fontsize',5);
    end;
    
    if plotnum==1
    %ylabel(strcat('Column:',int2str(con)));    
    ylabel(ytit,'Fontsize',5);
    end;
    set(gca,'Fontsize',5);
    %ylabel(strcat('Column ',cons));
    end
end

%text(-8,21.2,textdataDan(1).text,'units','centimeters');

%for j=1:nplots
    %dat=SerDan(j).SER;
%    TrapIntDatN2;
%end;
figlab=strcat('Averages of Time Series for column: ',cons);
%ccnplN;


%exname=strcat('c:\matlab6p1\work\',exdir,'tim');
%gcf=hf;
%print(gcf,'-djpeg','-r150',exname);