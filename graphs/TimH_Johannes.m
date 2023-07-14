%plots time height plots
fsize=18;
isave=0;
%load c:/matlabr12/work/bauru/casestudy/forcecons/diag/profnew+ccn_2-169;

minZ=0;
maxZ=22.58e3;  %19000;
maxtr=1.0;
time=Time;
tstart=12.67;
%dumprange=[50:104];
jmax=5; %max no. plots on one screen
a1=1;
a2=2; %values for subplot(ai,b,i)

iminovr=zeros([1 10]);
imaxovr=zeros([1 10]);
notsame=0; %flag to plot each plot individually in terms of colour scale


normcbar=0;

dumpint=300; %dump interval in seconds

i2d=0;
figlab='2d contour plot';

nplots=2;

ncont=25;
clab=0; %flag to label contours



clear h max timestxt pdat minC maxC



% for i=1:length(times)
%     te=num2str(times(i),3);
% 	timestxt(i,1:length(te))=te;
% end
% tit(1).tit='Max of Low Level Tracer';
% tit(2).tit='Max of Low Level Tracer';

logflag=0;


plotcase=1;

switch(plotcase)
    
	case 1
	fact=1;
    logflag=0;
    ncont=30;
    imaxovr(1:2)=1;
    maxcovOvr=500;
	tit(1).tit='Max Ice Number Concentration (/kg)';
	tit(2).tit='Max Ice Number Concentration (/kg)';
	
	case 2
	fact=1e3;
	tit(1).tit='Max Ice Mixing Ratio (g/kg)';
	tit(2).tit='Max Ice Mixing Ratio (g/kg)';
	
	case 3
	fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='Min Water Vapour Mixing Ratio (ppmv)';
	tit(2).tit=tit(1).tit;
	
	case 4
	fact=1;
	logflag=0;
	tit(1).tit='Max Upper Tracer Mixing Ratio (g/kg)';
	tit(2).tit='Max Upper Tracer Mixing Ratio (g/kg)';
    
    case 5
	fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='Low updraught 10th Percentile Water Vapour Mixing Ratio (ppmv)';
	tit(2).tit='High updraught 10th Percentile Water Vapour Mixing Ratio (ppmv)';
    
    ncont=30;
    
    case 6
	%fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='Low updraught 25th Percentile Ice Saturation Mixing Ratio (ppmv)';
	tit(2).tit='High updraught 25th Percentile Ice Saturation Mixing Ratio (ppmv)';
    
    case 7
	%fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='10th Percentile Ice Saturation Mixing Ratio (ppmv)';
	tit(2).tit=tit(1).tit;
    
    case 8
	%fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='95th Percentile potential temp (K)';
	tit(2).tit=tit(1).tit;
    
    case 9
    i2d=1;
	%fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='Ice Sat MR (ppmv)';
	tit(2).tit=tit(1).tit;
    figlab='LEM dump 85 ice sat MR';
    
    case 10
	logflag=1;
	tit(1).tit='Max Ice No. Conc (#/kg)';
	tit(2).tit=tit(1).tit;
    figlab='max ice NC time-height';
    
    case 11
	logflag=1;
    iminovr=1;
    mincovOvr=-3;
    imaxovr=1;
    maxcovOvr=log10(15);
	tit(1).tit='Low Updraught Case Max Ice Mixing Ratio (ppmv vapour equivalent)';
	tit(2).tit='Low Updraught Case Max Ice Mixing Ratio (ppmv vapour equivalent)'
    figlab='max ice MR time-height';
    
    case 12
	logflag=1;
    iminovr=1;
    mincovOvr=-8;
	tit(1).tit='Max Snow Mixing Ratio (g/kg)';
	tit(2).tit=tit(1).tit;
    figlab='max snow MR time-height';
    
    case 13
	logflag=1;
    iminovr=1;
    mincovOvr=-8;
	tit(1).tit='Max Graupel Mixing Ratio (g/kg)';
	tit(2).tit=tit(1).tit;
    figlab='max graupel MR time-height';
    
    case 14
	logflag=0;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-1.5;
	tit(1).tit='Max Vapour Deficit (ppmv)';
	tit(2).tit=tit(1).tit;
    figlab='Max Vapour Deficit time-height';
    
    case 15
	logflag=1;
    fact=1e6*28.97/18;
    iminovr=0;
    mincovOvr=-1.5;
	tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    figlab='Max Vapour time-height';
    
    case 16
	logflag=0;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-1.5;
	tit(1).tit='Max Vapour Deficit (ppmv)';
	tit(2).tit=tit(1).tit;
    figlab='Max Vapour Deficit time-height';
    %imaxovr=1;
    %maxcovOvr=-6.1
    
    case 17
	logflag=1;
    fact=1e6*28.97/18;
    %iminovr(1:2)=0;
    mincovOvr=-1.5;
	tit(1).tit='Low Updraught Case Max Water Vapour (ppmv)';
	tit(2).tit='High Updraught Case Max Water Vapour (ppmv)';
    figlab='Max Vapour time-height';
    %dumprange=[1:64];
    %nplots=1;
    %hrstartles=12.67;
    ncont=20;
    clab=1;
    
    case 18
	logflag=1;
    %fact=1e6*28.97/18;
%     iminovr=1;
%     mincovOvr=-1.5;
	tit(1).tit='Min Sat Vap MR (ppmv)';
	tit(2).tit=tit(1).tit;
    figlab='Min Sat MR';
    %imaxovr=1;
    %maxcovOvr=-6.1
    
    
    case 19
	logflag=0;
    fact=1e6*28.97/18;
%     iminovr=1;
%     mincovOvr=-1.5;
	tit(1).tit='Min Sat Vap MR (ppmv)';
	tit(2).tit='Min Vap MR (ppmv)';
    figlab='Min Sat + Vap MR';
    %imaxovr=1;
    %maxcovOvr=-6.1
    
    
    case 20
	logflag=1;
    fact=1e6*28.97/18;
%     iminovr=1;
%     mincovOvr=-1.5;
	tit(1).tit='Max Ice MR (g/kg)';
	tit(2).tit='Max Snow MR (g/kg)';
    figlab='Min Sat + Vap MR';
    iminovr=1;
    mincovOvr=-6.1
    %nplots=1;
    a1=2; %so plots as though are two graphs
    notsame=1;
    
    case 21
	logflag=1;
    fact=1e6*28.97/18;
     iminovr=1;
     mincovOvr=-5;
	tit(1).tit='Max Ice MR (g/kg)';
	tit(2).tit='Max Snow MR (g/kg)';
    figlab='Min Sat + Vap MR';
    %iminovr=1;
    %mincovOvr=-6.1
    %nplots=1;
    a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots=1;
    normcbar=0;
    
    case 22
	logflag=1;
    fact=1e6*28.97/18;
     iminovr=1;
     mincovOvr=-5;
	tit(1).tit='Max Ice MR (g/kg)';
	tit(2).tit='Max Snow MR (g/kg)';
    figlab='Min Sat + Vap MR';
    %iminovr=1;
    %mincovOvr=-6.1
    %nplots=1;
    a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots=1;
    normcbar=0;
    
    case 23
	logflag=0;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-1e-7;
	tit(1).tit='Mass Flux of Snow (kg/m^2/s)';
	tit(2).tit='Mass Flux of Snow (kg/m^2/s)';
    figlab='Snow Flux';
    imaxovr=1;
    maxcovOvr=0;
    %nplots=1;
    %a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots=1;
    %normcbar=0;
    
    case 24
	logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
	tit(1).tit='Fall Speed Flux of Ice (kg/m^2/s)';
	tit(2).tit='Fall Speed Flux of Ice (kg/m^2/s)';
    figlab='Fall Speed Ice Flux';
    imaxovr=1;
    maxcovOvr=1e-9;
    %nplots=1;
    %a1=2; %so plots as though are two graphs
    %notsame=1;
    
    %nplots=1;
    %normcbar=0;
    
    case 25
	logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
	tit(1).tit='Fall Speed Flux of Ice (kg/m^2/s)';
	tit(2).tit='Fall Speed Flux of Ice (kg/m^2/s)';
    figlab='Fall Speed Ice Flux';
    %imaxovr=1;
    maxcovOvr=30;
    %nplots=1;
    a1=2; %so plots as though are two graphs
    %notsame=1;
    
    nplots=1;
    %normcbar=0;
    
    dz=Grid.Z(izmin+1:izmax+1)-Grid.Z(izmin:izmax);
    
    case 26
	logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
	tit(1).tit='Net Flux of All Ice Species (kg/m^2/s)';
	tit(2).tit='Net Flux of All Ice Species (kg/m^2/s)';
    figlab='Net Ice Flux';
    %imaxovr=[0 1];
    maxcovOvr=0.5e-3;
    %nplots=1;
    a1=2; %so plots as though are two graphs
    notsame=1;
    
    ncont=16;
    
    %nplots=1;
    %normcbar=0;
    
    case 27
	logflag=0;
    fact=1e6*28.97/18;
    %iminovr=1;
    mincovOvr=-1e-7;
	tit(1).tit='Average Sublimation Rate of Ice (ppmv/s)';
	tit(2).tit='Average Sublimation Rate of Ice (ppmv/s)';
    figlab='Average Sublimation Rate of Ice';
    %imaxovr=1;
    maxcovOvr=6.5e-5;
    
    case 28
	logflag=1;
    fact=1e6*28.97/18;
    iminovr=1;
    mincovOvr=-2;
	tit(1).tit='Average Mixing Ratio of All Ice Species (ppmv vapour equivalent)';
	tit(2).tit='Average Mixing Ratio of All Ice Species (ppmv vapour equivalent)';
    figlab='Average Ice MR';
    imaxovr=1;
    maxcovOvr=0.6;
    ncont=17;
    
    

   
    
    
    
    
    
    
    
end
%timesTH=[hrstartles+((dumprange-1)*dumpint)/3600];
timesTH=Time(1:64)'+tstart;
timesTH=[timesTH timesTH(end)+(dumpint/3600.*(1:5))];
scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
hf=figure('position',posit,'name',figlab);

%nplots=length(prof);


    if nplots==1
        a=a1;
    else
        a=a2;
    end
    b=ceil(min(nplots,jmax)/2);

for  i=1:nplots
	%exname=strcat('c:/matlabr12/work/bauru/tracersjan2005/force+3_3th3qv/TracerTimH-',num2str(i));
	%xdat(i).x=time;
	%xdat(i).x=datenum(2004,2,24,hrstartles+floor(xdat(i).x/3600),60*(xdat(i).x/3600-floor(xdat(i).x/3600)),0);
	
	scrsz=get(0,'ScreenSize');
	posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
	
	iz=find(Grid.Z<=maxZ & Grid.Z>=minZ);
  
    izmax=iz(end);
    izmin=iz(1);
%     if length(iz)>=1
% 		iz=iz(1);
%     else
%         iz=length(Grid.Z);
%     end
	
	
	%pcolor(9+time./3600,Grid.Z(1:iz)./1000,maxLowTracer(i).prof(1:iz,47:80));hc=colorbar;%shading interp
	
    h(i).h=subplot(a,b,i);
    
    switch plotcase
        
     case 1
         switch i
         case 1
            zz(i).z=GridDan(1).Z(1:123)*150/123/1000;
            pdat(i).p=squeeze(max(c(:,:,:,1:64),[],1));
         case 2
            zz(i).z=GridDan(1).Z(1:izmax)/1000;
            pdat(i).p=squeeze(max(icenc(1).i(1:izmax,:,:),[],2))/1e6;
         end
     case 2
        pdat(i).p=prof(i).icemax(izmin:izmax,dumprange)*fact;
     case 3
        pdat(i).p=prof(i).vapmax(izmin:izmax,dumprange)*fact;
     case 4
        pdat(i).p=prof(i).uppTRmax(izmin:izmax,dumprange)*fact;
     case 5
        pdat(i).p=fact*pcents_vap(i).p(izmin:izmax,dumprange,2);
     case 6
        pdat(i).p=pcents_icemr(i).p(dumprange,izmin:izmax,4)';
     case 7
        pdat(i).p=pcents_icemr(i).p(dumprange,izmin:izmax,3)';
     case 8
       pdat(i).p=permute(pcents_potemp(i).p(izmin:izmax,dumprange,8),[2 1 3])';
     case 9
       pdat(i).p=satmr(i).s(izmin:izmax,:,85);
     case 10
       pdat(i).p=squeeze(max(icenc(i).i(izmin:izmax,:,dumprange),[],2));
     case 11
       pdat(i).p=squeeze(max(icemr(i).i(izmin:izmax,:,dumprange),[],2)) * fact/1000 ...
                    + squeeze(max(snowmr(i).i(izmin:izmax,:,dumprange),[],2)) * fact/1000 ;
     case 12
       pdat(i).p=squeeze(max(snowmr(i).i(izmin:izmax,:,dumprange),[],2));
     case 13
       pdat(i).p=squeeze(max(graupelmr(i).i(izmin:izmax,:,dumprange),[],2));
     case 14
       vv=fact*repmat(mean(vap(i).v(izmin:izmax,:,dumprange(1)),2),[1 length(dumprange)]);
       pdat(i).p=squeeze(min(fact*vap(i).v(izmin:izmax,:,dumprange),[],2))-vv;
     case 15
       pdat(i).p=fact*squeeze(max(vap(i).v(izmin:izmax,:,dumprange),[],2));
     case 16
       vv=fact*repmat(mean(vap(i).v(izmin:izmax,:,dumprange(1)),2),[1 length(dumprange)]);
       pdat(i).p=squeeze(pcents(i).p(dumprange,izmin:izmax,2))'-vv;
     case 17
       pdat(i).p=fact*squeeze(min(vap(i).v(izmin:izmax,:,dumprange),[],2));
     case 18
       pdat(i).p=squeeze(min(satmr(i).s(izmin:izmax,:,dumprange),[],2));
     case 19
         switch i
         case 1
            pdat(i).p=squeeze(min(satmr(1).s(izmin:izmax,:,dumprange),[],2));
         case 2
            pdat(i).p=fact*squeeze(min(vap(1).v(izmin:izmax,:,dumprange),[],2));
         end
    
     case 20
       switch i
         case 1
             pdat(i).p=squeeze(max(icemr(i).i(izmin:izmax,:,dumprange),[],2));
         case 2
            pdat(i).p=fact*squeeze(max(snowmr(1).i(izmin:izmax,:,dumprange),[],2));
         end
         
         case 21
            pdat(i).p=squeeze(max(V(i).v(izmin:izmax,:,dumprange),[],2));
            
        case 22
            %pdat(i).p=squeeze(max(Vsnow(i).v(izmin:izmax,:,dumprange),[],2));
            pdat(i).p=pcents_vsnow(i).p(izmin:izmax,dumprange,6);
            
         case 23
            %pdat(i).p=squeeze(max(Vsnow(i).v(izmin:izmax,:,dumprange),[],2));
            dgfind=findhead('ALL_WQ04',dgstrDan(1).dg)
            pdat(i).p=squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange));
            
         case 24
            pdat(i).p=squeeze(Falldiag(i).dg(izmin:izmax,6,dumprange));
            
         case 25
            dzz=repmat(dz,[1 length(dumprange)]);
            rho=repmat(GridDan(i).RHON(izmin:izmax),[1 length(dumprange)]);
            dgfind=findhead('ACC_A',dgstrDan(1).dg);
            A=squeeze(diag(i).dg(izmin:izmax,dgfind(1),1:length(dumprange)));
            az=find(A<0.001);
            A(az)=1;
    
            
            pdat(i).p=fact*squeeze(Falldiag(i).dg(izmin:izmax,6,dumprange))./dzz./rho./A*300;
            
          case 26
            pdat(i).p=squeeze(Fluxdiag(i).dg(izmin:izmax,6,dumprange) +Fluxdiag(i).dg(izmin:izmax,6+14,dumprange) - Falldiag(i).dg(izmin:izmax,6,dumprange))...
                + squeeze(Fluxdiag(i).dg(izmin:izmax,4,dumprange) +Fluxdiag(i).dg(izmin:izmax,4+14,dumprange) - Falldiag(i).dg(izmin:izmax,4,dumprange))...
                + squeeze(Fluxdiag(i).dg(izmin:izmax,5,dumprange) +Fluxdiag(i).dg(izmin:izmax,5+14,dumprange) - Falldiag(i).dg(izmin:izmax,5,dumprange));
            
        case 27
            pdat(i).p=fact*squeeze(pimlt(i).dg(izmin:izmax,1,dumprange) + psmlt(i).dg(izmin:izmax,1,dumprange) - pgmlt(i).dg(izmin:izmax,1,dumprange));
            
        case 28
            dgfind=findhead('ACC_A',dgstrDan(1).dg);
            A=squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange));
            az=find(A<0.001);
            A(az)=1;
                    
            dgfind=findhead('ALL_Q06',dgstrDan(1).dg);
            pdat(i).p=fact*squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange) )./A; 
            
            dgfind=findhead('ALL_Q04',dgstrDan(1).dg);
            pdat(i).p=pdat(i).p + fact*squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange) )./A; 
            
            dgfind=findhead('ALL_Q05',dgstrDan(1).dg);
            pdat(i).p=pdat(i).p + fact*squeeze(diag(i).dg(izmin:izmax,dgfind(1),dumprange) )./A; 
       
     end
    
        
    if logflag==1
        pdat(i).p=log10(pdat(i).p);
    end
    
    maxC(i)=max(max(pdat(i).p));
    minC(i)=min(min(pdat(i).p));
    

    
end





	maxCov=max(maxC);
	minCov=min(minC);






%need to make sure highest contour is higher than highest data value for colorbarf

if i2d==1
    timesTH=Grid.Y1/1000;
    xlabelstr='Horizontal Distance (km)';
else
    xlabelstr='Local Time (hrs)';
end


for i=1:nplots
    
    if notsame==1 
        maxCov=maxC(i);
	    minCov=minC(i);
    end
    
    if iminovr(i)==1
        minCov=mincovOvr;
    end

    if imaxovr(i)==1
        maxCov=maxcovOvr;
    end

    minc=minCov-0.1*abs(minCov);
    maxc=maxCov+abs(maxCov*0.1);
    
    conts=[minc:(maxc-minc)/ncont:maxc];

    ac=find(conts>0);
    if length(ac)>0 & ac(1)>1
        conts=[conts(1:ac-1) 0 conts(ac:end)];
    end

    
    half=abs((izmax-izmin)/2);
    pend=size(pdat(i).p,2);
    
    pdat(i).p(1:half,pend+1:pend+5)=minCov;
    pdat(i).p(half+1:end,pend+1:pend+5)=maxCov;
    
    h(i).h=subplot(a,b,i);
    
  
    
	%pcolor(timesTH,0.62+Grid.Z(izmin:izmax)./1000,pdat);
    
    [cbfA(i).c cbfB(i).c]=contourf(timesTH,zz(i).z,pdat(i).p,conts);
    if clab==1
        ch=cbfA(i).c;
        if logflag==1
            jc=1;
            while jc<size(cbfA(i).c,2)
                ch(1,jc)=str2num(num2str(10^(cbfA(i).c(1,jc)),'%2.0f'));
                jc=jc+cbfA(i).c(2,jc)+1;
            end
        end
        ch=round2(ch,0); %round to n decimal places
        clabel(ch,cbfB(i).c,'labelspacing',300);
    end
    %hc=colorbarf(cbfA,cbfB);
    shading flat
	
    %xti=set(h(i).h,'xticklabels',timestxt);
    
    set(h(i).h,'fontsize',fsize);
    xlabel(xlabelstr);
	ylabel('Height (km)');
    title(tit(i).tit);
    
    
    
    %caxis(h(i).h,[minCov maxCov*1.05]);
    if normcbar==1
        hc=colorbar( 'peer' , h(i).h );
	else
        hc=colorbarf(cbfA(i).c,cbfB(i).c); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!
    end
    
    if normcbar==0
    
      if logflag==1
        
        clear ctickstr
        
        ctick=get(hc,'yticklabel');
        ctickstr(1,1)=' ';
        for j=2:length(ctick)-1
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            nu=str2num(ctick{j});
            
                te=num2str(10^nu,3);       %'%2.2e');
   
            
            ctickstr(j,1:length(te))=te;
        end
        
        set(hc,'yticklabel',ctickstr);
        
%         add=str2num(ctick{end-1})/50;
%         
%         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
%         
%         for i=2:length(ctick)-1
%             cticknums(i)=str2num(ctick{i});
%         end
%         text(  ones( length(cticknums),1 )*1.05,cticknums+add,ctickstr, 'fontsize',fsize-6  );

      else
%         clear ctickstr
%         ctick=get(hc,'ytick');
%         for j=1:length(ctick)
%             %te=strcat('10^','{',num2str(ctick(j)),'}' );
%             te=num2str(ctick(j),'%2.2e');
%             ctickstr(j,1:length(te))=te;
%         end
%         
%         set(hc,'yticklabel','');
%         
%         add=ctick(end)/50;
%         
%         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
%         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
%     
%     
%         
%         set(hc,'yticklabel','');
%         
%         add=ctick(end)/50;
%         
%         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
%         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
%         
%         set(hc,'fontsize',fsize-6);
        
      end
    
    else

    set(hc,'fontsize',fsize-6);
    
    if logflag==1
        clear ctickstr
        ctick=get(hc,'ytick');
        for j=1:length(ctick)
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            te=num2str(10^ctick(j),3);
            ctickstr(j,1:length(te))=te;
        end
        
        set(hc,'yticklabel','');
        
        add=ctick(end)/50;
        
        set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
    
    
        
        set(hc,'yticklabel','');
        
        add=ctick(end)/50;
        
        set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
    
     else
        clear ctickstr
        ctick=get(hc,'ytick');
        for j=1:length(ctick)
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            te=num2str(ctick(j),'%2.2e');
            ctickstr(j,1:length(te))=te;
        end
     end
    
    
    
    set(hc,'fontsize',fsize-6);
    
    end
    
ylims=get(h(i).h,'ylim');   %re-scale to hide extra column put in to get the colorbars the same in both plots
axis(h(i).h,[timesTH(1) timesTH(end-5) ylims]);
    
end




% 
% tims=[9:2:23];
% ti=datenum(2004,2,24,tims,0,0);
% set(gca,'xtick',[ti]);
% datetick('x',15,'keepticks');

set(gcf,'paperpositionmode','auto');


if isave==1
     set(gcf,'paperpositionmode','auto');
	print(gcf,'-djpeg','-r350',exname);
    %print(gcf,'-dmeta',exname);
	%close(gcf);
end

