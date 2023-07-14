
scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.13];


f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
%  Ms=1/6*pi*rhoS*Ds.^3;
% 
% 
% a=find(times~=0);
% dt(2:length(times))=times(2:end)-times(1:end-1);
% dt(a(1))=dt(a(2));
% a=find(times==0);
% dt(a)=1; %avoid divide by zero
% dtt=repmat(dt,[150 1]);
% dtt2=repmat(dtt,[1 1 14]);
% dtt3=permute(dtt2,[1 3 2]);

xmin='';
xmax='';

if ~exist('justplot')
    justplot=0;
end

if justplot==0

clear labs xdat ydat
figname='Johannes graph';
gridon=1; %switch for grid default =1 - note probs with grid when resizing due to extra ticks.s
zmax='';
zmin=0;
lor=0;
noplot=0;
logflag=0;
lor=1;

xlab='x';
ylab='Height (km)';
     
     
    
graph=15;

switch graph
    case 15
    logflag=0;
    xlab='Temp (^oC)';
    ylab='Saturation over Ice'; %(cm^{-3})';
	
    

    figname=strcat('Sat over ice');
    
    v=4.3e-3;
    P=5e3;
    %P=5.54e3;
    T=-85:-70;
    
    PP=ones(size(T))*P;
    s=v./(SatVapPress(273.15-(70:85),'lem','ice',PP)/P * 18/28.79 *1e3);
    
    j=1;
    ydat(j).y = fliplr(s);
    xdat(j).x = T;
    labs(j).l='LEM';
    
    xmin=T(1);
    xmax=T(end);

    
    case 14
    logflag=0;
    xlab='Time UTC';
    ylab='Max Number Concentration (mg^{-1})'; %(cm^{-3})';
	
    

    figname=strcat('Max Number Concentration - All ice in LEM');
    
    j=1;
    ydat(j).y = maxNIlem + maxNSlem + maxNGlem
    xdat(j).x = time;
    labs(j).l='LEM';
    
    j=2;
    ydat(j).y = maxNImpc;
    xdat(j).x = time;
    labs(j).l = 'MPC';
    
    xmin=time(1);
    xmax=time(end);
    
   
  
 
    case 13
    logflag=0;
    xlab='Time UTC';
    ylab='Average Number Concentration (mg^{-1})'; %(cm^{-3})';
	
    
    tt=53;
    figname=strcat('Average Number Concentration - All ice in LEM');
    
    j=1;
    ydat(j).y = meanNIlem + meanNSlem + meanNGlem
    xdat(j).x = time;
    labs(j).l='LEM';
    
    j=2;
    ydat(j).y = meanNImpc;
    xdat(j).x = time;
    labs(j).l = 'MPC';
    
    xmin=time(1);
    xmax=time(end);
    
    
	case 12  
    logflag=0;
    xlab='Time UTC';
    ylab='Average Number Concentration of Ice Points (mg^{-1})'; %(cm^{-3})';
	
    
    figname=strcat('Average Number Concentration - All ice in LEM');
    
    j=1;
    ydat(j).y = meanNIlem2
    xdat(j).x = time;
    labs(j).l='LEM';
    
    j=2;
    ydat(j).y = meanNImpc2;
    xdat(j).x = time;
    labs(j).l = 'MPC';
    
    xmin=time(1);
    xmax=time(end);
    
    case 11
	logflag=1;
    
    xlab='Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='Beginning and end MAX and MIN profiles';
    
    
    
    
    xdat(1).x=pcents_icemr(2).p(50,2:end,1);
    ydat(1).y=Grid.Z(2:end)/1000;
    labs(1).l='High Updraught Beginning Ice Sat';
    

    j=2;
    xdat(j).x=pcents_icemr(2).p(end,2:end,1);
    ydat(j).y=Grid.Z(2:end)/1000;
    labs(j).l='High Updraught End Min Ice Sat';
    
    j=3;
    xdat(j).x=pcents_icemr(2).p(end,2:end,end);
    ydat(j).y=Grid.Z(2:end)/1000;
    labs(j).l='High Updraught End Max Ice Sat';
    
    j=4;
    xdat(j).x=f*mean(vap(2).v(:,:,50),2);
    ydat(j).y=Grid.Z(1:end)/1000;
    labs(j).l='High Updraught Beginning Vapour';
    

 
    
    case 10
          
        
    for i=1:2
        xdat((i-1)*2+1).x=f*max(TwoDDan(i).Q(:,2:end,1),[],2);
        ydat((i-1)*2+1).y=GridDan(i).Z/1000;
     
        xdat((i-1)*2+2).x=f*min(TwoDDan(i).Q(:,2:end,1),[],2);
        ydat((i-1)*2+2).y=GridDan(i).Z/1000;
        
    end
    
    T=tempLES(GridDan(1)); %K
    ei=SatVapPress(T,'goff','ice'); %Pa
    P=GridDan(1).PREFN; %Pa
    
    xdat(5).x=f*0.622*ei./(P-ei);
    ydat(5).y=GridDan(1).Z/1000
    
    T=tempLES(GridDan(2)); %K
    ei=SatVapPress(T,'goff','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    xdat(6).x=f*0.622*ei./(P-ei);
    ydat(6).y=GridDan(2).Z/1000;
    
    figname='Max/Min Vapour + Ice Sat';
        %ydat(1).y=(sumPosDep-sumPosUnDep)./sumPosDep;    
    
	labs(1).l='Max for Low Updraught Case';
   	labs(2).l='Min for Low Updraught Case';
	labs(3).l='Max for High Updraught Case';
	labs(4).l='Min for High Updraught Case';
    labs(5).l='Ice Sat Mixing Ratio Low Updraught';
    labs(6).l='Ice Sat Mixing Ratio High Updraught';
	
	xlab='Water Vapour Mixing Ratio (ppmv)';
    %xlab='Aerosol mass (kg)';
	ylab='Height (km)';
    
    
    logflag=1;
    lor=1;
    gridon=1;
    
    %set(gca,'xlim',[9e-9 2e-6]); do this after


    
    case 9
        
  
     xlab='Water Vapour Mixing Ratio (ppmv)';
     ylab='Height (km)';
     logflag=1;
     lor=1;
     
%      pcplot3(1,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'Low Updraught'},f);
%      
%      pcplot3(2,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'High Updraught'},f);
    
    pcs=[0 25 50 75 100];    

    [xdat,ydat,labs]=getpctdat(GridDan(2),f*TwoDDan(2).Q(:,2:end,1),pcs,' high updraught');    

    j=length(pcs)+1;
    
    
     xdat(j).x=f*pr(1).p(:,10);
     ydat(j).y=pr(1).p(:,1)/1000;
     labs(j).l='Bauru Sounding 17:15';
     
     xdat(j+1).x=f*pr(2).p(:,10);
     ydat(j+1).y=pr(2).p(:,1)/1000;
     labs(j+1).l='Campo Grande Sounding 9am';

     xdat(j+2).x=origVap;
     ydat(j+2).y=GridDan(1).Z/1000;
     labs(j+2).l='Original Sounding';

    case 8
    
    logflag=1;
    figname='Mean Percentiles Vapour Graph 50-104';
    xlab='Water Vapour Mixing Ratio (ppmv)';
    
    for i=50:size(vap(1).v,3)
        
        pcs=[0 25 50 75 100];
        
        for j=1:length(pcs)
            pcents(2).p(i,:,j)=prctile(f*vap(2).v(:,:,i)',pcs(j));
        end
    end
    
    for j=1:length(pcs)
        xdat(j).x=mean(pcents(2).p(:,:,j),1);
        ydat(j).y=GridDan(1).Z/1000;
    end
        
    
    lp=length(pcs);
    for j=1:lp
    
        if pcs(j)==0
            labstr='Min';
        elseif pcs(j)==100
            labstr='Max';
        else
            labstr=strcat( num2str(pcs(j)) , 'th percentile' );
        end
    labs(j).l=strcat(labstr,' for ',' high updraught case');
    end
    
    j=length(pcs)+1;
    xdat(j).x=origVap;
    ydat(j).y=GridDan(1).Z/1000;
    labs(j).l='Original Sounding';
    
     
    case 7
        
    logflag=1;

    xlab='Water Vapour Mixing Ratio (ppmv)';
    ylab='Height (km)';
    
%     for i=50:size(vap(2).v,3)
%         
%         pcs=[0 25 50 75 100];
%         
%         for j=1:length(pcs)
%             pc(i,:,j)=prctile(f*vap(2).v(:,:,i)',pcs(j));
%         end
%     end
    
    pcs=[0 25 50 75 100];

    for j=1:length(pcs)
        xdat(j).x=mean(pcents(1).p(50:104,:,j),1);
        ydat(j).y=620/1000+GridDan(1).Z/1000;
    end
        
    
    lp=length(pcs);
    for j=1:lp
    
        if pcs(j)==0
            labstr='Min';
        elseif pcs(j)==100
            labstr='Max';
        else
            labstr=strcat( num2str(pcs(j)) , 'th percentile' );
        end
    labs(j).l=strcat(labstr,' for ',' high updraught case');
    end
    
    j=length(pcs)+1;
    xdat(j).x=origorig;
    ydat(j).y=620/1000+GridDan(1).Z/1000;
    labs(j).l='Dump 1';
    
    j=length(pcs)+2;
    xdat(j).x=min(pcents(1).p(50:104,:,1),[],1); %minimum value over all times
    ydat(j).y=620/1000+GridDan(1).Z/1000;
    labs(j).l='Overall Min';
    
    j=length(pcs)+3;
    xdat(j).x=max(pcents(1).p(50:104,:,5),[],1); %max value over all times
    ydat(j).y=620/1000+GridDan(1).Z/1000;
    labs(j).l='Overall Max';
    

    
    
    case 6
        
    logflag=1;

    xlab='Water Vapour Mixing Ratio (ppmv)';
    
   
        xdat(1).x=f*mean(TwoDDan(1).Q(:,:,1),2);
        ydat(1).y=GridDan(1).Z/1000;
    
        labs(1).l='Start Water Vapour Profile';
        
        figname='Start Vapour Profile';
    
        
    case 5
     xlab='Water Vapour Mixing Ratio (ppmv)';
     ylab='Height (km)';
     logflag=1;
     lor=1;
     
%      pcplot3(1,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'Low Updraught'},f);
%      
%      pcplot3(2,1,[0 25 50 75 100],TwoDDan,GridDan,xlab,'',{'High Updraught'},f);
    
    %pcs=[0 25 50 75 100];    

    %[xdat,ydat,labs]=getpctdat(GridDan(2),f*TwoDDan(2).Q(:,2:end,1),pcs,' high updraught');    

    stats={'Corumba','CG','SP','Curitiba','Bauru'};

    for j=1:length(pro)
       
     xdat(j).x=f*pro(j).p(:,10,1);
     ydat(j).y=pro(j).p(:,1,1)/1000;
     labs(j).l=stats{j};
     
    end
    
    j=length(pro)+1;
%     xdat(j).x=f*mean(vap(1).v(:,:,50),2);
    xdat(j).x=origorig;
    ydat(j).y=620/1000 + GridDan(1).Z/1000; %Bauru 620m above msl
    labs(j).l='dump 1';
    
   
    
    case 4
        
    logflag=1;

    xlab='Water Vapour Mixing Ratio (ppmv)';
    
   
        
    figname='Dump 50 percentiles';
    
    pcs=[0 25 50 75 100];    

    [xdat,ydat,labs]=getpctdat(GridDan(2),f*vap(1).v(:,:,50),pcs,' high updraught');   

    
    
    
    
%     j=length(pro)+1;
%     xdat(j).x=f*mean(vap(1).v(:,:,50),2);
%     ydat(j).y=GridDan(1).Z/1000;
%     labs(j).l='dump 50';
    
    
case 3
    %run readsound, readSDLA, readallSAW
    
    logflag=1;
    
    xdat(1).x=sdla(7,:);
    ydat(1).y=sdla(2,:)/1000;
    labs(1).l='SDLA SF4 descent 18:57-21:48 LT';
    
   
    
    
    T=sdla(6,:)+273.15; %K
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=sdla(5,:)*100; %Pa
    
    xdat(2).x=f*0.622*ei./(P-ei);
    ydat(2).y=sdla(2,:)/1000;
    labs(2).l='Ice Sat Mixing Ratio SDLA';
    
    
%     T=tempLES(GridDan(2)); %K
%     ei=SatVapPress(T,'buck2','ice'); %Pa
%     P=GridDan(2).PREFN; %Pa
%     
%     xdat(3).x=f*0.622*ei./(P-ei);
%     ydat(3).y=GridDan(2).Z/1000;
%     labs(3).l='Ice Sat Mixing Ratio High Updraught';
    
  
    
%     
%     T=tempLES(GridDan(2)); %K
%     ei=SatVapPress(T,'goff','ice'); %Pa
%     P=GridDan(2).PREFN; %Pa
%     
%     xdat(4).x=f*0.622*ei./(P-ei);
    

    xdat(3).x=f*pr(1).p(:,10);
    ydat(3).y=pr(1).p(:,1)/1000;
    labs(3).l='DMI 17:15-19:08 LT';
    
    T=pr(1).p(:,3)+273.15; %K
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=pr(1).p(:,2)*100; %Pa
    
    xdat(4).x=f*0.622*ei./(P-ei);
    ydat(4).y=pr(1).p(:,1)/1000;
    labs(4).l='Ice Sat Mixing Ratio DMI';
    
    xdat(5).x=data(4).sss(7,:);
    ydat(5).y=data(4).sss(3,:)/1000;
    labs(5).l='SAW water vapour';
    
    %xdat(6).x=f*origorig;
    %ydat(6).y=620/1000 + GridDan(1).Z/1000; %Bauru 620m above msl
    %labs(6).l='LEM start profile';
    
    
    
    
    
	
    
    xlab='Water Vapour Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='SDLA profile + ice sat';
    
case 2
    logflag=1;
    
    xlab='Water Vapour Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='LEM anvil and outside profiles dump 85';
    
    x=3;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=1;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    x=20;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=2;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    x=40;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=3;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    x=60;
    xstr=num2str(GridDan(2).Y1(x));
    T=potemp(2).p(:,x,85)./(1e5./GridDan(2).PREFN).^0.286;
    
    ei=SatVapPress(T,'buck2','ice'); %Pa
    P=GridDan(2).PREFN; %Pa
    
    j=4;
    xdat(j).x=f*0.622*ei./(P-ei);
    ydat(j).y=GridDan(2).Z/1000;
    
    labs(j).l=strcat('Ice Sat Mixing Ratio x=',xstr,'  km');
    
    
    

case 1
	logflag=1;
    
    xlab='Water Vapour Mixing Ratio (ppmv)';
	ylab='Height (km)';
    
    figname='Beginning and end min profiles';
    
    
    xdat(1).x=min(vap(2).v(:,:,50),[],2);
    ydat(1).y=Grid.Z/1000;
    labs(1).l='High Updraught Beginning';
    
    j=2;
    xdat(j).x=min(vap(2).v(:,:,end),[],2);
    ydat(j).y=Grid.Z/1000;
    labs(j).l='High Updraught End';
    

end


end


if noplot==0
figure('name',figname,'Position',posit);
[h,ax,ax2]=plotXY4(xdat,ydat,labs,0,4,lor,logflag,xlab,ylab,[xmin xmax],[zmin zmax],1);
%function [H1,ax2]=plotXY3(xdat,ydat,labs,nmark,lwidth,leglor,logflag,xlab,ylab,ylims,zline) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers

% LEGEND(...,Pos) places the legend in the specified
%     location:
%         0 = Automatic "best" placement (least conflict with data)
%         1 = Upper right-hand corner (default)
%         2 = Upper left-hand corner
%         3 = Lower left-hand corner
%         4 = Lower right-hand corner
%        -1 = To the right of the plot

if gridon==1
    grid on;
end

end

