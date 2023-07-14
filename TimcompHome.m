%plots graphs from N direcs (NO diff graphs & six to a page)

%function [ha,hb]=timcomp(SerDan,direcDan,textdataDan,nplots,con);
clear rem;

profile=0;
pnames;

procon=0;

if ~exist('ovride');ovride=0;end
if ~exist('ovride2');ovride2=0;end
if ovride==0
    nplots=input('Enter number of plots :');
    con=input('Enter col no.: ');
end
con2=con;
scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

if ischar(con); profile=1; cons=con; end %flag for plotting profiles or timeseries (for profiles

clear minsx maxsx minsy maxsy sum i;
for iplot=1:nplots
    if ovride==1;i=iis(iplot);else i=iplot;end
    
    size(SerDan(i).SER(:,1));
    timsiz=ans(1);
    
    sgz=size(GridDan(i).Z,1);
    
    if profile==1
        ylabelDan='Height (km)';
    else
        xlabelDan='Time (hrs)';
    end
    
    if ovride==0;
        if profile==1
            xlabelDan=strcat('Column ',cons);
        else
            ylabelDan=strcat('Column ',cons);
        end
            titleDan{iplot}=direcDan(i).dir;
    end
        
           
   
    
 if ovride2==0 %if ==1 then put own pdata(i).p in
    if profile==1
        pdata2(iplot).p=GridDan(i).Z;
        %pdata(i).p=GridDan(i).OLQBAR(1:sgz,1)-Vbar(i).v(1:sgz);
        %pdata(i).p=meanW(i).w(1:sgz);
        if procon==0
            pdata(iplot).p=GridDan(i).VBAR(1:sgz,1);
        end
        if procon==1
        dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
        dgfind(1)
        pdata(iplot).p=firstdiag(i).DGAV(1:sgz,dgfind(1));
        end
        if procon==3
        dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
        dgfind(1)
        FirstGrid(i).OLQBAR(1,1)=0;
        pdata(iplot).p=diag(i).dg(1:sgz,dgfind(1))-firstdiag(i).DGAV(1:sgz,dgfind(1));
        end
        if procon==2
            dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
            ivap=dgfind(1);
            
            dgfind=findhead('ALL_TEMP',dgstrDan(i).dg);
            itemp=dgfind(1);
            
            Cp=1004.67*(1+0.84*diag(i).dg(1:sgz,ivap));
            Lv=(2.501-0.00237*(diag(i).dg(1:sgz,itemp)-273.16))*1e6;
            
            pdata(iplot).p= Cp.*diag(i).dg(1:sgz,itemp) + 9.81*GridDan(i).Z + Lv.*diag(i).dg(1:sgz,ivap);
        end
        
        if procon==4
            dgfind=findhead('ALL_WQ01',dgstrDan(i).dg);
            pdata(iplot).p=diagALL(i).dg(1:sgz,dgfind(1));%-firstdiag(i).DGAV(1:sgz,dgfind(1));
        end
        
        if procon==5
            dgfind=findhead('ALL_Q03',dgstrDan(i).dg);
            pdata(iplot).p=diagALL(i).dg(1:sgz,dgfind(1));%-firstdiag(i).DGAV(1:sgz,dgfind(1));
        end
        
    elseif con==999
        siser=size(SerDan(i).SER,1);
        totlat=0;
        pdata(iplot).p(1)=0;
        for iv=2:siser
            spacing=SerDan(i).SER(iv,1)-SerDan(i).SER(iv-1,1);
            totlat=totlat+0.5*spacing*(SerDan(i).SER(iv,2)+SerDan(i).SER(iv-1,2));
            vaplost=-SerDan(i).SER(iv,20)+SerDan(i).SER(1,20)+totlat/2.45e6;
            
            pdata(iplot).p(iv)=sum(SerDan(i).SER(iv,[21:25,56]))/vaplost;
            if vaplost<0.1;pdata(iplot).p(iv)=0;end
            
        end
        pdata(iplot).p=pdata(iplot).p';
        pdata2(iplot).p=SerDan(i).SER(:,1);
    else
        size(con);
        pdata(iplot).p=sum(  SerDan(i).SER(:,con(1:ans(2))),2  );
        pdata2(iplot).p=SerDan(i).SER(:,1);
    end
end %if ovride2==0
    
    
           
            
           
            
    
%     pdata(iplot).p=SerDan(i).SER(:,con(1));
%     size(con);
%     for icon=2:ans(2)
%         pdata(iplot).p=pdata(iplot).p+SerDan(i).SER(:,con(icon));
%     end
    
%     if con==0  %for addition plots
%         pdata(iplot).p=SerDan(i).SER(:,41)+SerDan(i).SER(:,42)+SerDan(i).SER(:,43);
%     else
%         pdata(iplot).p=SerDan(i).SER(1:timsiz,con);
%     end
    
    if profile==1
        minsx(iplot)=min(pdata(iplot).p);
        maxsx(iplot)=max(pdata(iplot).p);
        minsy(iplot)=min(pdata2(iplot).p);
        maxsy(iplot)=max(pdata2(iplot).p*1.1./1000);
	else
        minsx(iplot)=min(pdata2(iplot).p)./3600;
        maxsx(iplot)=max(pdata2(iplot).p)./3600;
        minsy(iplot)=min(pdata(iplot).p)/1.1.^sign(max(pdata(iplot).p));
        maxsy(iplot)=max(pdata(iplot).p)*1.1.^sign(max(pdata(iplot).p));
	end
    
end; %i=1:nplots
minx=min(minsx);
maxx=max(maxsx);
miny=min(minsy);
maxy=max(maxsy);

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
    
%     for idan=1:12
%         size(pdata(idan).p);
%         if ans(1)==296
%         else
%             break
%         end
%     end
    
    
    if profile==1
        figlab=strcat('Profiles of col: ',con);
    else
        cons=int2str(con);
        if con==108
            %figlab=strcat('Qcode: ',int2str(qcode),' v.time');
        elseif con>=60 & con<108
            figlab=strcat(pname(con).p,' v.time');
        else
            figlab=strcat('Column:',cons,' v.time');
        end
    end
    ha=figure('name',figlab,'Position',posit);
    set(ha,'paperpositionmode','auto');
    
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
        
        
        %if pro==1
          %  plot(pdata(j).p,SerDan(j).SER(:,1));
           % axis([miny maxy minx maxx]);
            % else
            if profile==1
                plot(pdata(j).p,pdata2(j).p./1000,'k');
            else
                plot(pdata2(j).p./3600,pdata(j).p,'k');
            end
            axis([minx maxx miny maxy]);
        %end
        
        
        title(titleDan{j},'fontsize',fsize);
          xlabel(xlabelDan,'fontsize',fsize);
          ylabel(ylabelDan,'fontsize',fsize);
          set(gca,'fontsize',fsize);
    end %j=(k-1)*jmax+1:min(jmax,nplots);
end %for k=1:fix((nplots-0.1)/nops)+1

if exist('textdataDan');text(-8,21.2,textdataDan(1).text,'units','centimeters');end
%if (profile==0 )%& nplots>1)
    
	for j=1:nplots
        clear dat;
        if profile==0
            dat(:,1)=pdata2(j).p;
            dat(:,2)=pdata(j).p;
        else
            dat(:,1)=pdata2(j).p;
            dat(:,2)=pdata(j).p.*GridDan(j).RHON;
        end
        
        
        con=2;
        TrapIntDatN2;
	end;
    if exist('textdataDan');
		if con2==108
            %figlab=strcat('Averages of Time Series for Qcode: ',int2str(qcode),',',textdataDan(1).text);
		elseif con2>=60
            figlab=strcat('Averages of Time Series for: ',pname(con).p,',',textdataDan(1).text);
		else
            figlab=strcat('Averages of Time Series for column: ',cons,',',textdataDan(1).text);
		end
    else
        figlab=strcat('Averages of Time Series for column: ',cons);
    end
    hb=figure('name',figlab,'Position',posit);
    set(hb,'paperpositionmode','auto');
	ccnplN;
    %end
con=con2;