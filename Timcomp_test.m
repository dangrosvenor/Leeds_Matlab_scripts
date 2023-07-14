%plots graphs from N direcs (NO diff graphs & six to a page)

lor=2; %1=right, 2=left, 0=automatic


clear rem llab H1 hs;

scrsz=get(0,'ScreenSize');
%posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
posit=[9 50 scrsz(3)/1.6 scrsz(4)/1.8]; %changed 06Mar08

if ~exist('ovride');ovride=0;end
if ~exist('iscript');iscript=0;end
if ~exist('noavplots');noavplots=0;end
if ~exist('fsize');fsize=12;end



profile=0; %or enter 'text' into col no.
pnames;



if ovride==0 
procon=1;
allonone=1;
end

if ovride==0
upheight=5.0; %range for profiles (km)
lowheight=0;
end

cdan(1).c=[1 0 0];  %red
cdan(3).c=[0 0.5 0.1]; %dark green
cdan(2).c=[0 0 1];  %blue
cdan(4).c=[0.5 0.5 0.5];  %grey
cdan(5).c=[0 0 0];
cdan(6).c=[0.7 0.8 0]; %yellowish

size(cdan);
scdan=ans(2);

pdan(1).p='-';
pdan(2).p=':';
%pdan(3).p='--';
%pdan(4).p='-.';

markers(1).m='none';
markers(2).m='+';
markers(3).m='o';
markers(4).m='*';
markers(5).m='.';
markers(6).m='x';
markers(7).m='s';
markers(8).m='d';
markers(9).m='^';
markers(10).m='<';
markers(11).m='>';
markers(12).m='p';
markers(13).m='h';
smark=length(markers);

linet=0;

spdan=length(pdan);


if ovride==0
    nplots=input('Enter number of plots :');
    con=input('Enter col no.: ');
end
con2=con;
con_save=con;
scrsz=get(0,'ScreenSize');

% if isstr(con)
% %     for idg=1:length(pname)
% %         if strcmp(pname(idg).p,con)==1
% %             con=idg;
% %             break
% %         end
% %     end
% 
% 
% end



dirs=1:nplots;
dirs2=nplots(2:end);
nplots=nplots(1);
if isempty(nplots);nplots=ndir;end


%if ischar(con); profile=1; cons=con; end %flag for plotting profiles or timeseries (for profiles

clear minsx maxsx minsy maxsy pdata sum i;
for ii=1:nplots
    
    if length(dirs2)>=1
        i=dirs2(ii);
    else
        i=ii;
    end
    
    if iscell(con_save)
     con=get_pname_col3(con_save,i,pname_all)
     con2=con;    
    end
    
    size(SerDan(i).SER(:,1));
    timsiz=ans(1);
    
    
    if upheight>0.01
        sgz=find(GridDan(i).Z./1000>upheight);
        sgz=sgz(1);
    else
        sgz=length(GridDan(i).Z);
    end
    
    if lowheight>0.01
        sgl=find(GridDan(i).Z./1000<lowheight);
        sgl=sgl(end);
    else
        sgl=1;
    end
    
    if ovride==0
        llab(ii).l=runName(i).nam;
    else
        llab(ii).l=tits(ii).tit;
    end
    
    
    if profile==1
        %pdata(i).p=GridDan(i).OLQBAR(1:sgz,1)-Vbar(i).v(1:sgz);
        %pdata(i).p=meanW(i).w(1:sgz);
        
        if iscript==1
            report2in;
        end
        
        otherdata(ii).dat=GridDan(i).Z(sgl:sgz)./1000;
        
        if procon==0
            pdata(ii).p=GridDan(i).OLQBAR(sgl:sgz,1);
        end
        if procon==1
            dgfind=findhead('ALL_WQ14',dgstrDan(i).dg);
            %dgfind(1)
            pdata(ii).p=TimeAvDan(i).DGAV(sgl:sgz,dgfind(1));
        end
        if procon==3
            dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
            dgfind(1)
            FirstGrid(i).OLQBAR(1,1)=0;
            pdata(ii).p=diag(i).dg(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
        end
        if procon==2
            dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
            ivap=dgfind(1);
            
            dgfind=findhead('ALL_TEMP',dgstrDan(i).dg);
            itemp=dgfind(1);
            
            Cp=1004.67*(1+0.84*diag(i).dg(sgl:sgz,ivap));
            Lv=(2.501-0.00237*(diag(i).dg(sgl:sgz,itemp)-273.16))*1e6;
            
            pdata(ii).p= Cp.*diag(i).dg(sgl:sgz,itemp) + 9.81*GridDan(i).Z + Lv.*diag(i).dg(sgl:sgz,ivap);
        end
        
        if procon==4
            dgfind=findhead('ALL_WQ01',dgstrDan(i).dg);
            pdata(ii).p=diag(i).dg(sgl:sgz,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
        end
        
        if procon==5
            dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
            pdata(ii).p=diag(i).dg(sgl:sgz,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
        end
        
         if procon==6
            dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
            pdata(ii).p=diag(i).dg(sgl:sgz,dgfind(1)).*GridDan(i).RHON(sgl:sgz);%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
%             dgfind=findhead('CLu_WQ05',dgstrDan(i).dg);
%             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1)).*GridDan(i).RHON(sgl:sgz);
%             
%             dgfind=findhead('CLu_WQ06',dgstrDan(i).dg);
%             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1)).*GridDan(i).RHON(sgl:sgz);
        end
        
        if procon==7
             dgfind=findhead('ALu_WQ04',dgstrDan(i).dg);
             pdata(ii).p=diag(i).dg(sgl:sgz,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
             
            dgfind=findhead('ALu_WQ05',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALu_WQ06',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALu_WQ01',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALu_WQ02',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALu_WQ03',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1));
        end
        
        if procon==8
           
             dgfind=findhead('ALL_Q04',dgstrDan(i).dg);
             pdata(ii).p=diag(i).dg(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
             
            dgfind=findhead('ALL_Q05',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALL_Q06',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q02',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q03',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+diag(i).dg(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
     
        end
        
        if procon==9
           
             dgfind=findhead('ALL_Q04',dgstrDan(i).dg);
             pdata(ii).p=firstdiag(i).DGAV(sgl:sgz,dgfind(1));
             
            dgfind=findhead('ALL_Q05',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALL_Q06',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q02',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q03',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+firstdiag(i).DGAV(sgl:sgz,dgfind(1));
     
        end
        
        if procon==10
           
             dgfind=findhead('ALL_Q04',dgstrDan(i).dg);
             pdata(ii).p=TimeAvDan(i).DGAV(sgl:sgz,dgfind(1));
             
            dgfind=findhead('ALL_Q05',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALL_Q06',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q02',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q03',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1));
     
        end
        
        if procon==11
           
             dgfind=findhead('ALL_Q04',dgstrDan(i).dg);
             pdata(ii).p=TimeAvDan(i).DGAV(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
             
            dgfind=findhead('ALL_Q05',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALL_Q06',dgstrDan(i).dg);
            pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q02',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
             dgfind=findhead('ALL_Q03',dgstrDan(i).dg);
             pdata(ii).p=pdata(ii).p+TimeAvDan(i).DGAV(sgl:sgz,dgfind(1))-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
     
        end
        
        if procon==12
           
             dgfind=findhead('ALL_Q01',dgstrDan(i).dg);
             %pdata(ii).p=GridDan(i).OLQBAR(sgl:sgz,1);
             pdata(ii).p=firstdiag(i).DGAV(sgl:sgz,dgfind(1));
         end
        
         if procon==14
             part='ACu';
             pdata(ii).p=zeros(1+sgz-sgz);
             for i14=2:6                           %2:6 = all condensate
                 nam=strcat(part,'_WQ0',int2str(i14));    
                 dgfind=findhead(nam,dgstrDan(i).dg);
                 %pdata(ii).p=GridDan(i).OLQBAR(sgl:sgz,1);
                 pdata(ii).p=pdata(ii).p+diag(ii).dg(sgl:sgz,dgfind(1));
             end
         end
         
    else
        size(con);
        pdata(ii).p=sum(  SerDan(i).SER(:,con(1:ans(2))),2  );
        otherdata(ii).dat=SerDan(i).SER(:,1)./3600;
    end
    
%     pdata(ii).p=SerDan(i).SER(:,con(1));
%     size(con);
%     for icon=2:ans(2)
%         pdata(ii).p=pdata(ii).p+SerDan(i).SER(:,con(icon));
%     end
    
%     if con==0  %for addition plots
%         pdata(ii).p=SerDan(i).SER(:,41)+SerDan(i).SER(:,42)+SerDan(i).SER(:,43);
%     else
%         pdata(ii).p=SerDan(i).SER(1:timsiz,con);
%     end
    
    if profile==1
        minsx(i)=min(pdata(ii).p);
        maxsx(i)=max(pdata(ii).p);
        minsy(i)=min(GridDan(i).Z(sgl:sgz))/1.1./1000;
        maxsy(i)=max(GridDan(i).Z(sgl:sgz))*1.1./1000;
	else
        minsx(i)=min(SerDan(i).SER(:,1))./3600;
        maxsx(i)=max(SerDan(i).SER(:,1))./3600;
        minsy(i)=min(pdata(ii).p)/1.1^sign(max(pdata(ii).p));
        maxsy(i)=max(pdata(ii).p)*1.1^sign(max(pdata(ii).p));
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

if minx>=maxx
    minx=-1;
    maxx=1;
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
    hfig=figure('name',figlab,'Position',posit);
    set(hfig,'paperpositionmode','auto');
    
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
        
        if allonone==0;
            subplot(a,b,plotnum);
            colour(j).c='b';
            patt(j).p='-';
        else
            hold on;
            if rem(j,scdan)==0
                colour(j).c=cdan(scdan).c;
            else
                colour(j).c=cdan(rem(j,scdan)).c;
            end
            if rem(j,spdan)==0
                patt(j).p=pdan(spdan).p;
            else
                patt(j).p=pdan(rem(j,spdan)).p;
            end
            if rem(j,smark)==0
                mark(j).m=markers(smark).m;
            else
                mark(j).m=markers(rem(j,smark)).m;
            end
            
            if rem(j,3)==0
                linet(j)=3;
            else
                linet(j)=rem(j,3);
            end
            
        end
        
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
                H1(j).h=plot(pdata(j).p,otherdata(j).dat);
            else
                H1(j).h=plot(otherdata(j).dat,pdata(j).p);
            end
            axis([minx maxx miny maxy]);
            
            if allonone==1
                markpnts=1:round(length(pdata(j).p)/20):length(pdata(j).p);
                if profile==1
                    H2(j).h=plot(pdata(j).p(markpnts),otherdata(j).dat(markpnts));
                    
                    
                else
                    H2(j).h=plot(otherdata(j).dat(markpnts),pdata(j).p(markpnts));
                   
                end
            
                set(H1(j).h,'LineStyle',patt(j).p);
                set(H1(j).h,'color',colour(j).c);
                set(H1(j).h,'linewidth',2);
                %set(H1(j).h,'linewidth',linet(j));
                set(H2(j).h,'LineStyle','none');
                set(H2(j).h,'color',colour(j).c);
                set(H2(j).h,'marker',mark(j).m);
                if strcmp(mark(j).m,'.')==1
                    set(H2(j).h,'markersize',20);
                end
                %set(H1(j).h,'marker',mark(j).m);
                
          
            end
        %end
        
        if allonone==0;title(llab(j).l);end
        
        if profile==1 & (allonone==0 | j==1)
            if ovride==0
                xlabel(cons,'fontsize',fsize);
            else
                xlabel(axlab,'fontsize',fsize);
            end
            ylabel('Height (km)','fontsize',fsize);
        elseif profile==0 & (allonone==0 | j==1)
            xlabel('Time (hrs)','fontsize',fsize);
        
            %ylabel('Precipitation (mm/hour)');
            if ovride==0
                if iscell(con_save)
                    ylabel(con_save,'fontsize',fsize);
                else                                        
                    ylabel(['Column ' cons],'fontsize',fsize);
                end
            else
                ylabel(axlab,'fontsize',fsize);
            end
        end
        
        set(hs(j).h,'fontsize',fsize);
        
    end %j=(k-1)*jmax+1:min(jmax,nplots);
end %for k=1:fix((nplots-0.1)/nops)+1
    
    if allonone==1
        if profile==1;
            [lh ha hb hc]=legend([H1.h],llab.l,lor);
        else
            [lh ha hb hc]=legend([H1.h],llab.l,lor);
        end
        set(lh,'Fontsize',fsize);
    end
    
 for k=1:fix((nplots-0.1)/nops)+1
     for j=(k-1)*jmax+1:min(jmax,nplots)
         set(ha(j*2+1),'marker',mark(j).m);
         set(ha(j*2+1),'color',colour(j).c);
         if strcmp(mark(j).m,'.')==1
            set(ha(j*2+1),'markersize',20);
         end
     end
 end
         
    
    
        
if exist('textdataDan');text(-8,21.2,textdataDan(1).text,'units','centimeters');end
%if (profile==0 )%& nplots>1)

if noavplots==0 & nplots~=1
	for j=1:nplots
        clear dat;
        if profile==0
            dat(:,1)=otherdata(j).dat;
            dat(:,2)=pdata(j).p;
        else
            dat(:,1)=otherdata(j).dat;
            dat(:,2)=pdata(j).p.*GridDan(j).RHON(sgl:sgz);
        end
        
        
        con=2;
        TrapIntDatN2;
	end;
    if exist('textdataDan');
		if con2==108
            %figlab=strcat('Averages of Time Series for Qcode: ',int2str(qcode),',',textdataDan(1).text);
		elseif con2>=60 & profile~=1
            figlab=strcat('Averages of Time Series for: ',pname(con).p,',',textdataDan(1).text);
        elseif profile==1
            figlab=strcat('Averages of Profiles of:',con2,',',textdataDan(1).text);
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

end