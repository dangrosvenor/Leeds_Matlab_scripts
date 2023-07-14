%plots graphs from N direcs (NO diff graphs & six to a page)


ownlab=1;
tit(1).t='Constant Fluxes';
tit(2).t='Sine wave flux variation';
tit(3).t='Double sine wave flux variation';

pnames;
nplots=input('Enter number of plots :');
con=input('Enter col no.: ');
con2=con;
scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];


exdir='field\cons\timcomp\';

clear minsx maxsx minsy maxsy pdata sum i;
for i=1:nplots
      
    size(SerDan(i).SER(:,1));
    timsiz=ans(1);
    
    size(con);
    pdata(i).p=sum(  SerDan(i).SER(:,con(1:ans(2))),2  );
    
%     pdata(i).p=SerDan(i).SER(:,con(1));
%     size(con);
%     for icon=2:ans(2)
%         pdata(i).p=pdata(i).p+SerDan(i).SER(:,con(icon));
%     end
    
%     if con==0  %for addition plots
%         pdata(i).p=SerDan(i).SER(:,41)+SerDan(i).SER(:,42)+SerDan(i).SER(:,43);
%     else
%         pdata(i).p=SerDan(i).SER(1:timsiz,con);
%     end
    
    
    minsx(i)=min(SerDan(i).SER(:,1));
    maxsx(i)=max(SerDan(i).SER(:,1));
    minsy(i)=min(pdata(i).p);
    maxsy(i)=max(pdata(i).p);
end;
minx=min(minsx)./3600;
maxx=max(maxsx)./3600;
miny=min(minsy);
maxy=max(maxsy)*1.1;

if miny>=maxy
    miny=-1;
    maxy=1;
end



exist posit;
a=ans;
if a==0;
    posit=[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2.2];
end;

jmax=12;
nops=12; %no plots per screen
for k=1:fix((nplots-0.1)/nops)+1
    
    for idan=1:12
        size(pdata(idan).p);
        if ans(1)==296
        else
            break
        end
    end
    
    
        cons=int2str(con);
        if con==108
            figlab=strcat('Qcode: ',int2str(qcode),' v.time');
        elseif con>=60 & con<108
            figlab=strcat(pname(con).p,' v.time');
        else
            figlab=strcat('Column:',cons,' v.time');
        end
    
        
    
    hff.f=figure('name',figlab,'Position',posit);
    
    
    for j=(k-1)*jmax+1:min(jmax,nplots);
    plotnum=rem(j,jmax);
    if plotnum==0;
        plotnum=jmax;
    end
    
    if nplots==1
        a=1;
    else
        a=2;
    end
    b=ceil(min(nplots,jmax)/2);
    
    a=3;
    b=1;
    
    
    subplot(a,b,plotnum);
    
%     jmax=nops;
%     if nplots-k*nops<0 
%         jmax=rem(nplots,nops)+(k-1)*nops;
%     end;
%     
%     for j=(k-1)*nops+1:jmax*k;
%         plotnum=rem(j,nops);
%         if plotnum==0;
%             plotnum=nops;
%         end
    
    %subplot(6,2,plotnum);
    
    hs(j).h=newplot;
    
    
    
    plot(SerDan(j).SER(:,1)./3600,pdata(j).p);
    axis([minx maxx miny maxy]);
    if ownlab==0
        title(direcDan(j).dir);
    else
        title(tit(j).t);
    end
    
    if j==nplots
        xlabel('Time (hrs)');
    end
    
    ylabel('Max Height of Ice (m)');
    
    %ylabel(strcat('Column ',cons));
    end
end

text(-8,21.2,textdataDan(1).text,'units','centimeters');

for j=1:nplots
    clear dat;
    dat(:,1)=SerDan(j).SER(:,1);
    dat(:,2)=pdata(j).p;
    con=2;
    TrapIntDatN2;
end;
if con2==108
    figlab=strcat('Averages of Time Series for Qcode: ',int2str(qcode),',',textdataDan(1).text);
elseif con2>=60
    figlab=strcat('Averages of Time Series for: ',pname(con).p,',',textdataDan(1).text);
else
    figlab=strcat('Averages of Time Series for column: ',cons,',',textdataDan(1).text);
end


            
			
            exname=strcat('c:\matlabR12\work\',exdir,'timser',textdataDan(1).text,'-',int2str(k),'col-',cons,'.jpg');
			gcf=hff.f;
			set(gcf,'paperpositionmode','auto');
			print(gcf,'-djpeg','-r350',exname);
			close(gcf);

       
ccnplN;
gcf=hff.f;
set(gcf,'paperpositionmode','auto');

exname=strcat('c:\matlabR12\work\',exdir,'timAV',textdataDan(1).text,'-',int2str(k),'col-',cons,'.jpg');
			
			
			print(gcf,'-djpeg','-r150',exname);
			close(gcf);
