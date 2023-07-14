flag=21;

clear pdata pdata2 profstr
calc=1;
itim=1;
isave=0;
maxheight=-1;
iA=0;
factor=1;
imultDEN=0;
gridon=0;
idiff=0;

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

ovride=1;
ovride2=1;

nplots=3;
profile=0;
fsize=14;

style{1}='-';
style{2}='-.';
style{3}=':';
style{4}='--';

colour(1).c='k';
colour(2).c=[0 0.5 0.1];
colour(3).c=[0 0 1];
colour(4).c=[0.5 0.5 0.5];
width=[1 4 4 4]; %widths don't work for anything apart from '-' lines - work with exporting as JPEGS but not on screen!
%style{i}='k-';

saveDir='bauru\casestudy\graphsFrance\';

switch flag
    case 1
    titleDan{1}='13th Feb, 18:00';
    titleDan{2}='24th Feb, 20:15';
    titleDan{3}='3rd Mar, 12:00';
    ylabelDan='Cloud Top (m)';
    cols=[13 17 17];
    iis=[2 8 9];
    for i=1:nplots
        pdata(i).p=max(SerDan(iis(i)).SER(:,[13 17]),[],2);
        pdata2(i).p=SerDan(iis(i)).SER(:,1);
    end
    fname='CtopNorm';
    calc=0;
    isave=1;
    
    
    case 2
    titleDan{1}='13th Feb, 18:00';
    titleDan{2}='24th Feb, 20:15';
    titleDan{3}='3rd Mar, 12:00';
    xlabelDan='Time (hrs)';
    ylabelDan='Precipitaion (mm)';
    cols=[56];
    iis=[2 8 9];
    itim=0;
    fname='PrecipNorm';
    isave=1;
    
    case 3
    titleDan{1}='13th Feb, 18:00';
    titleDan{2}='24th Feb, 20:15';
    titleDan{3}='3rd Mar, 12:00';
    xlabelDan='Time (hrs)';
    ylabelDan='Max Updraught (m/s)';
    cols=[4];
    iis=[2 8 9];
    itim=0;
    fname='UpdraughtNorm';
    isave=1;
    
    case 4
    titleDan{1}='13th Feb, 18:00 + 3degC';
    titleDan{2}='24th Feb, 20:15 + 3degC';
    titleDan{3}='3rd Mar, 12:00 + 3degC';
    ylabelDan='Cloud Top (m)';
    cols=[13 17 17];
    iis=[1 5 3];
    for i=1:nplots
        pdata(i).p=max(SerDan(iis(i)).SER(:,[13 17]),[],2);
        pdata2(i).p=SerDan(iis(i)).SER(:,1);
    end
    calc=0;
    isave=1;
    fname='CloudTopsAlt';
    
    case 5
    titleDan{1}='13th Feb, 18:00 + 3degC';
    titleDan{2}='24th Feb, 20:15 + 3degC';
    titleDan{3}='3rd Mar, 12:00 + 3degC';
    ylabelDan='Precipitation (mm)';
    cols=[56];
    iis=[1 5 3];
    calc=1;
    isave=1;
    itim=0;
    fname='PrecipAlt';
    
    case 6
    titleDan{1}='13th Feb, 18:00 + 3degC';
    titleDan{2}='24th Feb, 20:15 + 3degC';
    titleDan{3}='3rd Mar, 12:00 + 3degC';
    ylabelDan='Updraught (m/s)';
    cols=[4];
    iis=[1 5 3];
    calc=1;
    isave=1;
    itim=0;
    fname='Updraught';
    
    case 7
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    ylabelDan='Updraught (m/s)';
    cols=[4];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    fname='Updraught24';
    
    case 8
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    ylabelDan='Precipitation (mm)';
    cols=[56];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    fname='Precip24';
    
    case 9
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    ylabelDan='Cloud Top (m)';
    cols=[13 17 17];
    iis=[8 6 5 7];
    for i=1:nplots
        pdata(i).p=max(SerDan(iis(i)).SER(:,[13 17]),[],2);
        pdata2(i).p=SerDan(iis(i)).SER(:,1);
    end
    calc=0;
    isave=1;
    itim=0;
    fname='Ctop24';
    
    case 10
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    ylabelDan='Downdraught (m/s)';
    cols=[5];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    fname='Downdraught24';
    
    case 11
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    ylabelDan='Precipitation Rate(mm/hr)';
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    fname='PrecipRate24';
    
    case 12
    nplots=4;
    ttitleDan{1}='0^{\circ}C';
    titleDan{2}='2^{\circ}C';
    titleDan{3}='3^{\circ}C';
    titleDan{4}='4^{\circ}C';
    xlabelDan='Water Flux(kg/m^2/s)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='ProfileWaterFlux24';
    profstr{1}='ALL_WQ01';
    profstr{2}='ALL_WQ02';
    profstr{3}='ALL_WQ03';
    profstr{4}='ALL_WQ04';
    profstr{5}='ALL_WQ05';
    profstr{6}='ALL_WQ06';
    
    maxheight=15; %height to scale vert axis to in km
    imultDEN=1;
    
    case 13
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    xlabelDan='Average Updraught(m/s)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='ProfileAluW24';
    profstr{1}='ALu_W';
    
    %maxheight=21.5; %height to scale vert axis to in km
    
    case 14
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    xlabelDan='Microphysical Source of Potential Temp (K/day)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='DTH';
    profstr{1}='ALL_DTH';
    maxheight=15;
    factor=3600*24;
    imultDEN=0;
    %maxheight=21.5; %height to scale vert axis to in km
    
    case 15
    nplots=4;
    titleDan{1}='0degC';
    titleDan{2}='2degC';
    titleDan{3}='3degC';
    titleDan{4}='4degC';
    xlabelDan='Average Updraught(m/s)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=0;
    itim=0;
    profile=1;
    fname='ProfileCluW24A';
    profstr{1}='CLu_W';
    iA=1;
    
    case 16
    nplots=4;
    titleDan{1}='0^{\circ}C';
    titleDan{2}='2^{\circ}C';
    titleDan{3}='3^{\circ}C';
    titleDan{4}='4^{\circ}C';
    xlabelDan='Microphysical Source of Potential Temp (K/day)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='LW24';
    profstr{1}='ALL_LW';
    maxheight=15;
    
    case 17
    nplots=4;
    titleDan{1}='0^{\circ}C';
    titleDan{2}='2^{\circ}C';
    titleDan{3}='3^{\circ}C';
    titleDan{4}='4^{\circ}C';
    xlabelDan='Heating due to SW+LW Radiation (K/day)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='LW+SW24';
    profstr{1}='ALL_SW';
    profstr{2}='ALL_LW';
    maxheight=15;
    
    case 18
    nplots=4;
    titleDan{1}='0^{\circ}C';
    titleDan{2}='2^{\circ}C';
    titleDan{3}='3^{\circ}C';
    titleDan{4}='4^{\circ}C';
    xlabelDan='Increase in Temperature ({\circ}C)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='IncreaseTemp24';
    profstr{1}='ALL_TEMP';
    maxheight=15;
    gridon=1;
    idiff=1;
    
    case 19
    nplots=4;
    titleDan{1}='0^{\circ}C';
    titleDan{2}='2^{\circ}C';
    titleDan{3}='3^{\circ}C';
    titleDan{4}='4^{\circ}C';
    xlabelDan='Change in Water Vapour (kg/kg)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='IncreaseVap24';
    profstr{1}='ALL_Q01';
    maxheight=15;
    gridon=1;
    idiff=1;
    
    case 20
    nplots=4;
    titleDan{1}='0^{\circ}C';
    titleDan{2}='2^{\circ}C';
    titleDan{3}='3^{\circ}C';
    titleDan{4}='4^{\circ}C';
    ylabelDan='Total Ice Hydormeteors (kg/m^2)';
    cols=[25];% 24 25];
    iis=[8 6 5 7];
    calc=1;
    isave=0;
    itim=0;
    fname='ISG24';
    
    case 21
    nplots=4;
    titleDan{1}='0^{\circ}C';
    titleDan{2}='2^{\circ}C';
    titleDan{3}='3^{\circ}C';
    titleDan{4}='4^{\circ}C';
    xlabelDan='Average Ice (kg/m^3)';
    ovride2=0;
    cols=[26];
    iis=[8 6 5 7];
    calc=1;
    isave=1;
    itim=0;
    profile=1;
    fname='IceProf24';
    profstr{1}='ALL_Q06';
    maxheight=15;
    gridon=1;
    idiff=0;
    imultDEN=1;
    
    
end

if calc==1
for i=1:nplots
    
    pdata(i).p=0;
    
    lc=length(cols);
    ilc=min(i,lc);
    if profile==0
        pdata(i).p=sum(SerDan(iis(i)).SER(:,cols(:,lc)),2);
        pdata2(i).p=SerDan(iis(i)).SER(:,1);
    elseif idiff==0
        sgz=size(diagALL(iis(i)).dg,1);
        for isum=1:length(profstr)
            dgfind=findhead(profstr{isum},dgstrDan(iis(i)).dg);
            dgsum(isum)=dgfind(1);
        
            if iA==1
                findA=findhead(strcat(profstr{isum}(1:4),'A'),dgstrDan(iis(i)).dg);
                divA(isum).d=diagALL(iis(i)).dg(1:sgz,findA(1));
            else
                divA(isum).d=1;
            end
            divA(isum).d=divA(isum).d/factor;
            if imultDEN==1;divA(isum).d=divA(isum).d./GridDan(iis(i)).RHON(1:sgz);end
            fprintf(1,'%d %s %d %s\n',dgfind(1),dgstrDan(iis(i)).dg{dgfind(1)},findA(1),dgstrDan(iis(i)).dg{findA(1)});
            
            iiA=find(divA(isum).d==0);
            divA(isum).d(iiA)=1;
            
            
            pdata(i).p=pdata(i).p+diagALL(iis(i)).dg(1:sgz,dgsum(isum))./divA(isum).d;
            pdata2(i).p=GridDan(iis(i)).Z(1:sgz);
        end %isum
        
    else %idiff==1
        for i=1:nplots
        sgz=size(diagALL(iis(i)).dg,1);
           
                dgfind=findhead(profstr{1},dgstrDan(iis(i)).dg);
                dgsum(1)=dgfind(1);
         
            
            pdata(i).p=diag(iis(i)).dg(1:sgz,dgsum(1))-firstdiag(iis(i)).DGAV(1:sgz,dgsum(1));
            pdata2(i).p=GridDan(iis(i)).Z(1:sgz);
        end
        
        
    end
end   
end 

if itim==1
    timcompHome;
else
    ha=figure('Position',posit);
    set(ha,'paperpositionmode','auto');
    for i=1:nplots
        if profile==0    
            plot(pdata2(i).p./3600,pdata(i).p,style{i},'linewidth',width(i),'color',colour(i).c);
            xlabelDan='Time (hrs)';
        else
            plot(pdata(i).p,pdata2(i).p./1000,style{i},'linewidth',width(i),'color',colour(i).c);
            ylabelDan='Height (km)';
        end
        
        if maxheight>0; set(gca,'ylim',[0 maxheight]); end
        xlabel(xlabelDan,'fontsize',fsize);
        ylabel(ylabelDan,'fontsize',fsize);
        set(gca,'fontsize',fsize);
        hold on;
        
        %title(titleDan{i});
        leg=legend(titleDan,0);
        set(leg,'fontsize',fsize);
    end
end

if gridon==1; grid; end

if isave==1
    
    print(ha,'-djpeg','-r100',strcat(saveDir,fname,'.jpg'));
    %close(ha);
end

ovride=0;
ovride=0;