%plots multiple line graphs with dual axis option
clear cdan AX exax;

timdiv=1; %factor for changing time to minutes=60 or hours=3600 seconds=1

cdan(1).c=[1 0 0];  %red
cdan(2).c=[0 0.5 0.1]; %dark green
cdan(4).c=[0 0 1];  %blue
cdan(5).c=[0.5 0.5 0.5];  %grey
cdan(3).c=[0 0 0];
cdan(6).c=[0.7 0.8 0]; %yellowish

size(cdan);
scdan=ans(2);

pdan(1).p='-';
pdan(2).p=':';



titsbaurufull;

tit(1).tit='Droplet conc. = 60 /cm^3';
tit(2).tit='Droplet conc. = 360 /cm^3';
tit(3).tit='Droplet conc. = 960 /cm^3';

clear conm llab;

nplots=input('Enter number of plots :');

group=10; %switch for showing different sets of results
dualax=1; %switch =1 to use dual axis graphs
expo=0; %switch to export file as jpeg 
exdir='bauru\homogen\processrates\ISGmrsJust3-2'; %filepath to export to
fontnam='times-roman';

lorr(1)=1;  %switches for left or right axis 0=left 1=right - put 1st right hand axis line in lorr(1)
lorr(2)=0;        
lorr(3)=0;
lorr(4)=0;
lorr(5)=0;
lorr(6)=0;

laxlab='Total Process Rate (kg/m^2/s)';
raxlab='Total Mixing Ratio (kg/m^2)';




if expo==1
    fsize=10;  %fontsize - choose 7 if you want to read it, 3 if exporting to file
else
    fsize=10;
end

if group==10
    conm(1).c=88;
    conm(2).c=90;
    
    conm(3).c=89;
    conm(4).c=91;
    %conm(5).c=92;
    %conm(6).c=91;
    %conm(7).c=93;
    
    %llab(1).l='Ice No conc' %44
    llab(2).l='Ice Deposition'; %90
    llab(3).l='Primary Nucleation '; %89
    llab(1).l='Hallet Mossop'; %88
    %llab(4).l='Contact Nucleation'; %92
    llab(4).l='Ice Accreting Water'; %91
    %llab(5).l='Freezing'; %93
    
    
    %lorr(1)=1;
    
    %laxlab='';
    %raxlab='Max No. Conc. (/kg)';
    dualax=1;
end

    
if group==9
    
    
    
    conm(1).c=44;
    conm(2).c=45;
    conm(3).c=46;
    
    llab(1).l='Max No. Conc. of Ice (/kg)';
    llab(2).l='Max No. Conc. of Graupel (/kg)';
    llab(3).l='Max No. Conc. of Snow (/kg)';
    
    laxlab='Max No. Conc. (/kg)';
    raxlab='Max No. Conc. (/kg)';
end
if group==8
    conm(1).c=69;
    conm(2).c=79;
    
    llab(1).l='IACR-G';
    llab(2).l='PGFR';
 
end
if group==7
    conm(1).c=23;
    conm(2).c=6;
    
    llab(1).l='Rain';
    llab(2).l='W flux';

elseif group==6
%showing total sink of liq

conm(1).c=21;
conm(2).c=[91 72 77 64 68 60 90];

llab(1).l='Liq';
llab(2).l='Sinks of liq.';

elseif group==1

%showing depletion of liquid water - depletion for low CCN due to rain acc. water & from ISG acc. water for high CCN - more water left for high CCN
conm(1).c=21;           
conm(2).c=[91 72 77];
conm(3).c=64;
%conm(4).c=[68 60 90];
conm(4).c=62;
%conm(6).c=108;

llab(1).l='Liq. W.';
llab(2).l='I+S+G Acc. W.';
llab(3).l='R Acc. W.';
%llab(4).l='I+S+G Dep.';
llab(4).l='Auto-con';
%llab(6).l='dQcond/dt';



elseif group==2
%showing different types of graupel formation - from rain for low CCN & from water and snow for high
conm(1).c=21;        
conm(2).c=[73 78]; %does sum of cols 73 & 78
conm(3).c=69;
conm(4).c=72;

llab(1).l='Liq W';
llab(2).l='G Auto+G Acc S';
llab(3).l='I Acc.R-G';
llab(4).l='G Acc.W';



elseif group==3
%showing different rain production modes - from auto and acc. low CCN and from graupel high CCN
conm(1).c=23;
conm(2).c=[61 65 66]; 
conm(3).c=62;
%conm(5).c=63;
conm(4).c=64;

llab(1).l='Total Rain MR';
llab(2).l='Melting I+S+G';
llab(3).l='Auto-con';
%llab(5).l='G Shedding';
llab(4).l='R Acc. W';


elseif group==4
%showing amounts of ice, snow and graupel

laxlab='Tot.MR(snow&graupel) (kg/m^2)';
raxlab='Tot.MR (ice) (kg/m^2)';


conm(1).c=22;
conm(2).c=24; 
conm(3).c=25;

llab(1).l='Ice';
llab(2).l='Snow';
llab(3).l='Graupel';
    
elseif group==5
%showing different rain production modes - from auto and acc. low CCN and from graupel high CCN

raxlab='Total Mixing Ratio (kg/m^2)';

conm(1).c=23;
conm(2).c=62;
conm(3).c=63;
conm(4).c=85

llab(1).l='Rain';
llab(2).l='Auto-con';
llab(3).l='G Shedding';
llab(4).l='Rain Evp.';


end





scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
hfig=figure('name','Combination plot','position',posit,'name',strcat('Combination Plot of group:',int2str(group)));




size(conm);
sconm=ans(2);
pnames;
%nplots=12;

%con=input('Enter col no.: ');
%fact=input('Please enter scale factor: ');


scrsz=get(0,'ScreenSize');

clear minsx maxsx minsy maxsy H1;
maxovall=0;
maxovallR=0;  %for right hand axis
clear pser;

for i=1:nplots
    size(conm);
    sicm=ans(2);
    maxxx(i)=max(SerDan(i).SER(:,1));
    
    
    
    for mc=1:sicm
        size(conm(mc).c);
        sisum=ans(2);
        for isum=1:sisum   %do sum for lines with e.g. [62 63 64]
            if isum==1
                pser(i).ser(:,mc)=SerDan(i).SER(:,conm(mc).c(isum));
            else
                pser(i).ser(:,mc)=pser(i).ser(:,mc)+SerDan(i).SER(:,conm(mc).c(isum));
            end
        end
        
        if lorr(mc)==0
            if max(pser(i).ser(:,mc))>maxovall
                maxovall=max(pser(i).ser(:,mc));
            end
        else
            if max(pser(i).ser(:,mc))>maxovallR
                maxovallR=max(pser(i).ser(:,mc));  %max for right hand axis
            end
        end
    end
end;
maxy=maxovall*1.1;
maxyR=maxovallR*1.1;
maxx=max(maxxx)/timdiv;


exist posit;
a=ans;
if a==0;
    posit=[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2.2];
end;

jmax=12;
for k=1:1
    %cons=int2str(con);
    %figure('name',strcat('Column:',cons,' v.time'),'Position',posit);
    
    
    
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
    
    %hold on;
    subplot(a,b,plotnum);
    hold on;
    
    hs(j).h=newplot;
    
    times=SerDan(j).SER(:,1)./timdiv;
    
    cdanc=0;
    pdanc=1;
    
    
    if dualax==1
        [AX H1(2).h H1(1).h]=plotyy(times,pser(j).ser(:,2),times,pser(j).ser(:,1));
        cmstart=3;
    else
        cmstart=1;
    end
    
    for mc=cmstart:sicm
        cdanc=cdanc+1;
        pdanc=pdanc+1;
        if cdanc==scdan+1;
            cdanc=1;
        end
        if pdanc==3
            pdanc=1;
        end
        
        size(SerDan(j).SER(:,1));
        stim=ans(1);
        size(pser(j).ser(:,mc));
        spser=ans(1);
        
        if lorr(mc)==0    %if is to be normally plotted on left
            H1(mc).h = plot(times,pser(j).ser(1:min(stim,spser),mc));
        else
            axes(AX(2));  %if is to be plotted on right
            hold on;
            H1(mc).h = plot(times,pser(j).ser(1:min(stim,spser),mc));
            axes(AX(1));  %reset back to left as AX may not exist for lorr=0
            hold on;
        end
        
        set(H1(mc).h,'LineStyle',pdan(pdanc).p);
        set(H1(mc).h,'color',cdan(cdanc).c);
        set(H1(1).h,'linewidth',1.7);
        if dualax==1
        set(H1(2).h,'linewidth',1.7);
         set(H1(3).h,'linewidth',1.7);
     end
    end
    
    if j==12   %add legend only on last plot
        
    end
    lh=legend([H1.h],llab.l,2);
    set(lh,'Fontsize',8);
    set(lh,'fontname',fontnam);
   
    exist AX;
    if ans==1
        exax=1;
        set(AX,'xlim',[0 maxx]);
        set(AX(1),'ylim',[0 maxy]);
        
        
        set(get(AX(1),'Ylabel'),'String',laxlab,'fontsize',fsize,'fontname',fontnam);
        %if max(lorr==1)
            set(get(AX(2),'Ylabel'),'String',raxlab,'fontsize',fsize,'fontname',fontnam);
            %end
        set(AX,'Fontsize',fsize);
        set(AX,'fontname',fontnam);
        
        if j==11 | j==12
            
            
            
        end
        
        ex=floor(log10(maxy)); %exponent
        nu=maxy/10^ex;         %number
        if nu<=2
            spa=0.5*10^ex;
        elseif nu<=5
            spa=1*10^ex;
        else
            spa=2*10^ex;
        end
       
        %set tick marks cuz Matlab makes them go weird!?
        yts=[0 spa spa*2 spa*3 spa*4];
        set(AX(1),'ytick',yts);
    
        %maxyR=8e5;
    %if max(lorr)==1
        set(AX(2),'ylim',[0 maxyR]);
        ex=floor(log10(maxyR)); %exponent
        nu=maxyR/10^ex;         %number
        if nu<=2
            spa2=0.5*10^ex;
        elseif nu<=5
            spa2=1*10^ex;
        else
            spa2=2*10^ex;
        end

        yts2=[0 spa2 spa2*2 spa2*3 spa2*4];
        set(AX(2),'ytick',yts2);
        
        %end
        
        
    else
        set(hs(j).h,'xlim',[0 maxx]);
        set(hs(j).h,'ylim',[0 maxy]);
        set(hs(j).h,'fontsize',fsize);
        set(hs(j).h,'fontname',fontnam);
        exax=0;
    end
    %axis([minx maxx miny maxy]);
    
    if j==nplots 
    xlabel('Time (mins)','Fontsize',fsize+1,'fontname',fontnam);
    end
    %xlabel('Time (s)','Fontsize',fsize+1,'fontname',fontnam,'position',[maxx/1.1 -maxy/5]);
    title(tit(j).tit,'Fontsize',fsize+1,'fontname',fontnam);
    
    
    end
end



if expo==1
    exname=strcat('c:\matlabR12\work\',exdir); 
    set(hfig,'paperpositionmode','auto');
    print(hfig,'-djpeg','-r150',exname);
    close(hfig);
end