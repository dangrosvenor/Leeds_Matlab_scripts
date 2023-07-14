%plots time height plots
fsize=18;
isave=0;
%load c:/matlabr12/work/bauru/casestudy/forcecons/diag/profnew+ccn_2-169;

minZ=0;
maxZ=25000;
maxtr=1.0;
time=[50*300:300:104*300];
hrstartles=9;
dumprange=[50:104];

dumpint=300; %dump interval in seconds



clear h max timestxt

timesTH=[9+((dumprange-2)*dumpint)/3600];

% for i=1:length(times)
%     te=num2str(times(i),3);
% 	timestxt(i,1:length(te))=te;
% end
% tit(1).tit='Max of Low Level Tracer';
% tit(2).tit='Max of Low Level Tracer';

logflag=0;


plotcase=5;

switch(plotcase)
    
	case 1
	fact=1;
	tit(1).tit='Max Ice Number Concentration (/kg)';
	tit(2).tit='Max Ice Number Concentration (/kg)';
	
	case 2
	fact=1e3;
	tit(1).tit='Max Ice Mixing Ratio (g/kg)';
	tit(2).tit='Max Ice Mixing Ratio (g/kg)';
	
	case 3
	fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='Max Water Vapour Mixing Ratio (ppmv)';
	tit(2).tit=tit(1).tit;
	
	case 4
	fact=1;
	logflag=0;
	tit(1).tit='Max Upper Tracer Mixing Ratio (g/kg)';
	tit(2).tit='Max Upper Tracer Mixing Ratio (g/kg)';
    
    case 5
	fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='Min Water Vapour Mixing Ratio (ppmv)';
	tit(2).tit=tit(1).tit;

end


scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
hf=figure('position',posit);

%nplots=length(prof);
nplots=2;

if nplots==1
        a=1;
    else
        a=2;
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
    izmin=iz(1)
%     if length(iz)>=1
% 		iz=iz(1);
%     else
%         iz=length(Grid.Z);
%     end
	
	
	%pcolor(9+time./3600,Grid.Z(1:iz)./1000,maxLowTracer(i).prof(1:iz,47:80));hc=colorbar;%shading interp
	
    h(i).h=subplot(a,b,i);
    
    switch plotcase
        
     case 1
        pdat=prof(i).iceNCmax(izmin:izmax,dumprange)*fact;
     case 2
        pdat=prof(i).icemax(izmin:izmax,dumprange)*fact;
     case 3
        pdat=prof(i).vapmax(izmin:izmax,dumprange)*fact;
     case 4
        pdat=prof(i).uppTRmax(izmin:izmax,dumprange)*fact;
     case 5
        pdat=pcents(i).p(dumprange,izmin:izmax,1)';
        
    end
        
    if logflag==1
        pdat=log10(pdat);
    end
    
    maxC(i)=max(max(pdat));
    
	%pcolor(timesTH,0.62+Grid.Z(izmin:izmax)./1000,pdat);
    
    [cbfA cbfB]=contourf(timesTH,0.62+Grid.Z(izmin:izmax)./1000,pdat,23);
    %hc=colorbarf(cbfA,cbfB);
    shading flat
	
    %xti=set(h(i).h,'xticklabels',timestxt);
    
    set(h(i).h,'fontsize',fsize);
    xlabel('Local Time (hrs)');
	ylabel('Height (km)');
    title(tit(i).tit);
    
    
    
end




maxCov=max(maxC);

% clear ctickstr
% for i=1:nplots
%     
%     set(gcf,'currentaxes',h(i).h);
%     %
%     
%     caxis(h(i).h,[0 maxCov*1.05]);
%     %hc=colorbar( 'peer' , h(i).h );    
%     hc=colorbarf(cbfA,cbfB);
%     
%     
%     
%     if logflag==1
%         ctick=get(hc,'ytick');
%         for j=1:length(ctick)
%             %te=strcat('10^','{',num2str(ctick(j)),'}' );
%             te=num2str(10^ctick(j),3);
%             ctickstr(j,1:length(te))=te;
%         end
%         
%         set(hc,'yticklabel','');
%         
%         add=ctick(end)/50;
%         
%         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
%         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
% 	end
%     
%     set(hc,'fontsize',fsize-6);
%     
% end




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

