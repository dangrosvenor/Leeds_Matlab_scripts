%imports several diagnostic files, plots 2d, saves as jpeg and closes down figure

titlecb='';
sorthead=0; %flag to switch off sorting of headers to save time
same=0; %flag for giving all plots the same colour scale (=1)
field=0; %flag for baurufield to add directories to filename
fsize=15; %fontsize
plotflag=0;
maxZovr=23e3; %set to zero for top of domain
minZovr=15e3;
noplot=1; %flag to stop figures appearing (=1)
output=0; %flag for output to file (=1)
contoura=[0 0]; %array to tell which plot to use contours on
samesize=0; %flag to plot all 2d plots on the same size domain for comparing runs with different domain sizes
samefile=0; %flag to allow plotting from a stored array

zlim=[1 0];
lims=[]; %user defined color limits for each plot lims(j).l(1) and lims(j).l(2)

limcol=0;
limcols=[7350 7650];

itext=0;
ititle=1; %flag to display titles
ititovr=1; %flag to override auto titles with manually typed ones

inoread=0; %flag so that goes through jj loop but doesn't read in data

logcol=0; %flag to plot log of colour values
fvap=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

useozmin=1; %flag to use value below as min value on colour scale for ozone plots
ozmin=5.3e-8;

nb=20;

anim=0;


iread=1; %flag to say whether to read in files or just plot the one that's already in memory
streamflag=0; %flag for putting streamlines of wind on
vectorf=1;

%exdirA='field\etasound\30.01-16\etavars\Radar\'; %backslash at end
exdirA='c:\matlabR12\work\bauru/tracersjan2005/aerosols/Nenes_HeightVar/droplets'; % / at end !
exdirA='c:/cygwin/home/login/runs/aeroruns/aeroqadv2_4_deplete/results/droplets/'
exdirA='c:\matlabR12\work\bauru/tracersjan2005/force+3_3th3qv/totalMR+iceSat_3';
exdirA='c:\matlabR12\work\bauru/tracersjan2005/sdlaVap/TotalMR';
%exdirA='c:\matlabR12\work\bauru/tracersjan2005/force+3_3th3qv/snowNC/';
exdirA='g:runs/dmi1715_5ppmv/results/icenc+cloud/';

if strcmp(exdirA(end),'/')==0; exdirA(end+1)='/'; end

clear b hh time bintime dav bav v rain w rem contour bin dav qf maxims minims;

fnmin=input('Enter the first diagnostic file number, with directories in direcDan(i).dir: ');
fnmax=input('Enter the maximum diagnostic file number, with directories in direcDan(i).dir: ');
morecon=input('Enter col no.: ');
cono=morecon;
nplots=input('Enter the no. of directories to import :');

if isstr(cono)==1 %if are using cono as a string to identify a process then set cono to an unused number
    constr=cono;
    cono=0;
    morecon=0;
else
    constr='';
end

dirs=nplots(2:end);
nfiles=nplots(1);
nplots=nplots(1);

moreplots=size(morecon,2);

if moreplots>1
    same=0;
    nfiles=1;
else
    nfiles=nplots;
    contours=contoura(1);
end

if iread==0 | fnmin==0
    nfiles=0;
    fnmin=2;
    fnmax=2;
    output=0;
    if length(dirs)>0
        for i=1:length(dirs)
            titles(i).tit=direcDan(dirs(i)).dir;
        end
    else
        for i=1:nplots
            titles(i).tit=direcDan(i).dir;
        end
    

    end
end

if ititovr==1
        %titles(1).tit='Ice Saturation Mixing Ratio (ppmv)';
        %titles(2).tit='High Updraught Case';
        %titles(2).tit='Total Condensate (g/kg)';
        
        titles(1).tit='';
        titles(2).tit='';
end

if samefile==1; nfiles=1; moreplots=2; morecon=ones([1 10]); end %if want to plot from values stored in an array


if cono==266
    qf=input('Enter Q field no. : ');
    qfstr=int2str(qf);
else
    qfstr='';
end

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
titsbaurufull;

if sorthead==1
sortheaderdgs;
%ndth=strmatch('ALL_DTH',dgstr); %column no. in TimeAv.DGAV for latent heat profiles due to microphysics
end

if inoread==1; nfiles=0; end



for jj=fnmin:fnmax  %for1
    fprintf(1,'jj=%d ',jj);
    if jj<10;
        fn=strcat('00',int2str(jj));
    elseif jj<100
        fn=strcat('0',int2str(jj));
    else
        fn=int2str(jj);
    end

    for i=1:nfiles
        
        if length(dirs)>=1 %if want to import only certain files then enter no. files and then file numbers to be imported
            idir=dirs(i); %e.g. enter [3 5 6 7] to enter files 5,6&7
        else
            idir=i;
        end
        
        %titles(i).tit=direcDan(idir).dir;
       
        if strcmp(direcDan(idir).dir(2:2),':')==1
            FileName=strcat(direcDan(idir).dir,'\RUN0001.DG0',fn);
        else
            FileName=strcat('c:\cygwin\home\user\',direcDan(idir).dir,'\RUN0001.DG00',fn);
        end
        
        DIAG_3d;
        GridDan(idir)=Grid;
        TwoDDan(idir)=TwoD;
        SerDan(idir).SER=SER;
        textdataDan(idir).text=textdat;
        %TimeAvDan(i)=TimeAv;
        try 
            TimeAvDan(idir).surav=TimeAv.surav;
        catch
        end   
        TimeAvDan(idir).RNDGS=TimeAv.RNDGS;
        TimeAvDan(idir).DGAV=TimeAv.DGAV;
        
        if jj==fnmin
            bin(idir).b(1:nb)=0;
            dav(idir).d(1:size(Grid.Y1,1))=0;
        end
    end;
  
sb=(size(Grid.Y1)-2)/nb; %size of bin for divergences
%time(jj-fnmin+1)=SER(end,1); %store times of each file

clear maxims rem;   

jmax=6;
kmax=fix((nplots-0.1)/6)+1;

    if moreplots>1; nplots=moreplots; end
    
    if nplots==1
        a=1;
    else
        a=2;
    end
    b=ceil(min(nplots,jmax)/2);
   

   
for k=1:kmax     %for2 no. figures
    
    
    
    if noplot==0
        %hh((jj-2)*kmax+k).h=figure('name',strcat('2d plots with same scales of column:',int2str(cono),' :',qfstr),'Position',posit);
        hff=figure('name',strcat('2d plots with same scales of column:',int2str(cono),' :',qfstr),'Position',posit);
    end
    
     if (nplots-k*6<0)
         jmax=rem(nplots,6)+(k-1)*6;
     end;

    jmax=nplots;

for j=(k-1)*jmax+1:min(jmax*k,nplots); %for3 number of subplots per page
 
    if moreplots>1
        cono=morecon(j);
        jc=1;
        contours=contoura(j);
    elseif length(dirs)>=1
        jc=dirs(j);
        cono=morecon;
    else
        cono=morecon;
        jc=j;
    end
    
    plotnum=rem(j,jmax);
    if plotnum==0;
        plotnum=jmax;
    end
   
    if (noplot==0)
        hs((jj-1)*nplots+j).h=subplot(a,b,plotnum); %subplotting
    end
    
    [ZZ XX]=meshgrid(GridDan(jc).Z,GridDan(jc).Y1);
    ZZ=ZZ';
    XX=XX';
    xxmax(j)=max(GridDan(jc).Y1);
    zzmax(j)=max(GridDan(jc).Z);
    
    
    
    %for j=1:nplots
    if (cono<=9 & cono>0)
        
        pdat=TwoDDan(jc).Q(:,:,cono);
        
        if cono==1;
            titlecb='Vapour MR (kg/kg)';
            pdat(:,1)=0;
        end;
        
        if cono==2;
            titlecb='Liquid MR (kg/kg)';
            pdat(:,1)=0;
        end;
        
        if cono==3;
            titlecb='Rain Water MR (kg/kg)';
        end;
        if cono==4;
            titlecb='Snow MR (kg/kg)';
        end;
        if cono==5;
            titlecb='Graupel MR (kg/kg)';
        end;
        if cono==6;
            titlecb='Ice MR (kg/kg)';
        end;
        if cono==7;
            titlecb='Ice NC (#/kg)';
        end;
        if cono==8;
            titlecb='Graupel NC (/kg)';
        end;
        if cono==9;
            titlecb='Snow NC (#/kg)';
        end;
        
    
    
        
    end %if cono<=9
    
    if strcmp(constr,'Total Aerosol')==1
        pdat=sum(TwoDDan(jc).Q(:,:,16:end),3)/1e6; %sum of all aerosol particles (cm^-3)
        titlecb=constr;
    end
        
    if cono==555
        for i555=1:size(TwoDDan(jc).TH1,1)
            pdat(i555,:)=TwoDDan(jc).TH1(i555,:)+GridDan(jc).THREF(i555);
        end
        contours=1;
        titlecb='';
    end;    
    
    
        
    
    if cono==1212
        pdat=TwoDDan(jc).Q(:,:,12);
        titlecb='';
    end;
    
    if cono==112
%          for ih=1:size(Grid.Z,1)
%              pdat(ih,:)=TwoDDan(jc).V(ih,:)-Grid.VBAR(ih);
%         end
        pdat=TwoDDan(jc).V;
        titlecb='V (m/s)';
    end;
    
    
    if cono==1121
        DYY=Grid.Y1(2)-Grid.Y1(1);
        for ih=1:size(Grid.Z,1)
            pdat(ih,:)=gradient(TwoDDan(jc).V(ih,:)-Grid.VBAR(ih),DYY);
        end
        titlecb='Divergence (s^-1)';
    end;
    
    if cono==1122
        DYY=Grid.Y1(2)-Grid.Y1(1);
        for ih=1:size(Grid.Z,1)
            pdat(ih,:)=gradient(TwoDDan(jc).V(ih,:)-GridDan(jc).VBAR(ih),DYY);
        end
        for ih=1:size(Grid.Z,1)
            divp(j).d(jj-fnmin+1,ih)=mean(pdat(ih,:)); %profile of horizontally averaged divergence at each time
        end
        
        dav1=mean(pdat(1:15,:)); %average over Grid.Z(15)=1.05km
        divs(j).d(:,:,jj-fnmin+1)=pdat; %save all data to save going through all files in future
        
        
        %dav(j).d=dav(j).d+dav1;  %sum divergences for all points in Y
        
        
%         
%         for ih=1:nb
%             bav(j).b(ih)=mean(dav1((ih-1)*sb+1:ih*sb)); %average over bin
%             bin(j).b(ih)=bin(j).b(ih)+bav(j).b(ih); %total value for this bin
%         end
%         bintime(j).b(jj-fnmin+1)=sum(bav(j).b); %sum over all bins for timeseries
%             plotflag=1;
%             %binall=binall+bin(ih).b./(fnmax-fnmin+1);
%         titlecb='Averaged & Binned (s^-1)';

    end
    
    if cono==1123
        %v(j).dat(:,:,jj)=TwoDDan(jc).V;
        %precip(j).dat(:,jj)=TwoDDan(jc).Q(1,:,3);
        for ih=1:size(Grid.Y1,1)
            rain(j).dat(ih,jj)=mean(TwoDDan(jc).Q(:,ih,3),1);
        end
        %w(j).dat(:,:,jj)=TwoDDan(jc).W;
        
    end
    
    if strcmp(constr,'SW rad')==1
        titlecb='SW rad';
        pdat=TwoDDan(jc).swrad;
    end;
    
    if strcmp(constr,'STH rad')==1
        titlecb='STH rad';
        pdat=TwoDDan(jc).sthrad;
    end;
    
    if strcmp(constr,'TWD rad')==1
        titlecb='TWD rad';
        pdat=TwoDDan(jc).twdrad;
    end;
    
    if cono==26;
        titlecb='Total LRSGI MR (g/kg)';
        pdat=1000*TwoDDan(jc).Q(:,:,2);
        for kk=3:6;
            pdat=pdat+1000*TwoDDan(jc).Q(:,:,kk); %in g/kg
        end;
        
    end;
    
    if cono==266

        if qf==133
            y1=-16.2;
            y0=-35;
            x1=exp(y1);
            if useozmin==0
                x0=min(min(TwoDDan(jc).Q(:,:,end)));
            else
                x0=ozmin;
            end
            mq=(y1-y0)/(x1-x0);
            kq=y1-(mq*x1);
            pdat=TwoDDan(jc).Q(:,:,14);
            aq=find(pdat<exp(y1));
            pdat(aq)=mq*pdat(aq) + kq;
            bq=find(TwoDDan(jc).Q(:,:,end)>=exp(y1));
            pdat(bq)=log(pdat(bq));

              % pdat=log(TwoDDan(jc).Q(:,:,13));
          elseif qf==1
              pdat=fvap*TwoDDan(jc).Q(:,:,qf); %vapour in ppmv
          else
            pdat=TwoDDan(jc).Q(:,:,qf);
        end
        titlecb=strcat('MR of Q field: ',num2str(qf),' (/kg)');
    end
    
    if cono==23
        noplot=1;
        output=0;
%         if jj==fnmin
%             diag(j).dg=zeros(size(TimeAvDan(j).DGAV));
%         end
        
        times(jj)=SerDan(j).SER(end,1);
       
       dgstore=[]; %store of all the indexes for partitioned properties  
       dgfind=findhead('ALL_A',dgstrDan(j).dg);
        for idg=1:27
            temp=dgstrDan(j).dg{dgfind(1)+idg-1}; %look through strings from ALL_A
            part=temp(1:3); %=partition names
            if idg~=6       %6=weird partitition that is just a blank
                dgfind2=findhead(part,dgstrDan(j).dg); %look for all DGAV no.s with that partitition
                dgstore=[dgstore dgfind2];
                dgfindArea=findhead( strcat(part,'_A') , dgstrDan(j).dg); %find area no. for this partition
                areacol=TimeAvDan(j).DGAV(:,dgfindArea(1));
                zeroareas=find(areacol==0);
                areacol(zeroareas)=1; %set to one to avoid divide by zero (DGAV values will be zero for these anyway)
                area=repmat(areacol,size(TimeAvDan(j).DGAV(1,dgfind2(2:end)))); %make 2-d matrix of replicas of area
                diag(j).dg(:,dgfind2(2:end),jj)=TimeAvDan(j).DGAV(:,dgfind2(2:end))./area;
                diag(j).dg(:,dgfind2(1),jj)=TimeAvDan(j).DGAV(:,dgfind2(1)); %don't divide by areas for areas themselves
            end     %'2:end' so don't include xxx_A
            
            dg=[1:size(TimeAvDan(j).DGAV,2)];
            dg(dgstore)=[]; %get rid of those already done above
            diag(j).dg(:,dg,jj)=TimeAvDan(j).DGAV(:,dg);  %puts them into diag
            
        end
            
      end
      
      if cono==232
        noplot=0;
        output=1;
        if jj==fnmin
            diag(j).dg=zeros(size(TimeAvDan(j).DGAV));
        end
        
        dgfind=findhead('ALL_A',dgstrDan(j).dg);
        for idg=1:27
            temp=dgstrDan(j).dg{dgfind(1)+idg-1}; %look through strings from ALL_A
            part=temp(1:3); %=partition names
            if idg~=6       %6=weird partitition that is just a blank
                dgfind2=findhead(part,dgstrDan(j).dg); %look for all DGAV no.s with that partitition
                dgfindArea=findhead( strcat(part,'_A') , dgstrDan(j).dg); %find area no. for this partition
                areacol=TimeAvDan(j).DGAV(:,dgfindArea(1));
                zeroareas=find(areacol==0);
                areacol(zeroareas)=1; %set to one to avoid divide by zero (DGAV values will be zero for these anyway)
                area=repmat(areacol,size(TimeAvDan(j).DGAV(1,dgfind2(2:end)))); %make 2-d matrix of replicas of area
                diag(j).dg(:,dgfind2(2:end),jj)=diag(j).dg(:,dgfind2(2:end))+TimeAvDan(j).DGAV(:,dgfind2(2:end))./area;
            end     %2:end so don't include xxx_A
            
            diag(j).dg(:,dgfind2(1),jj)=diag(j).dg(:,dgfind2(1))+TimeAvDan(j).DGAV(:,dgfind2(1)); %for areas themselves don't divide by areas
            
        end
        
        titlecb='Total LRSGI MR (/kg)';
        pdat=TwoDDan(jc).Q(:,:,2);
        for kk=3:6;
            pdat=pdat+TwoDDan(jc).Q(:,:,kk);
        end;
            
      end
    
    if cono==233
            dgfind=findhead('CLu_WQ01',dgstrDan(jc).dg);
            Clu_WQxx(j).prof(:,jj-fnmin+1,1)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('CLu_WQ02',dgstrDan(jc).dg);
            Clu_WQxx(j).prof(:,jj-fnmin+1,2)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('CLu_WQ03',dgstrDan(jc).dg);
            Clu_WQxx(j).prof(:,jj-fnmin+1,3)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('CLu_WQ04',dgstrDan(jc).dg);
            Clu_WQxx(j).prof(:,jj-fnmin+1,4)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('CLu_WQ05',dgstrDan(jc).dg);
            Clu_WQxx(j).prof(:,jj-fnmin+1,5)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('CLu_WQ06',dgstrDan(jc).dg);
            Clu_WQxx(j).prof(:,jj-fnmin+1,6)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            
            dgfind=findhead('ALu_WQ01',dgstrDan(jc).dg);
            Alu_WQxx(j).prof(:,jj-fnmin+1,1)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALu_WQ02',dgstrDan(jc).dg);
            Alu_WQxx(j).prof(:,jj-fnmin+1,2)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALu_WQ03',dgstrDan(jc).dg);
            Alu_WQxx(j).prof(:,jj-fnmin+1,3)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALu_WQ04',dgstrDan(jc).dg);
            Alu_WQxx(j).prof(:,jj-fnmin+1,4)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALu_WQ05',dgstrDan(jc).dg);
            Alu_WQxx(j).prof(:,jj-fnmin+1,5)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            dgfind=findhead('ALu_WQ06',dgstrDan(jc).dg);
            Alu_WQxx(j).prof(:,jj-fnmin+1,6)=TimeAvDan(jc).DGAV(:,dgfind(1));%-firstdiag(i).DGAV(sgl:sgz,dgfind(1));
            
            
    end 
    
    if cono==2332
            
%             col_str='ACu_W';
%             
%             area_str=strcat(colstr(1:3),'_A');
%             dgfind=findhead(area_str,dgstrDan(jc).dg);
%             area=TimeAvDan(jc).DGAV(:,dgfind(1));
%             areazeroes=find(area==0);
%             area(areazeroes)=1;
%             
%             dgfind=findhead(col_str,dgstrDan(jc).dg);
%             ACu_W(j).prof(:,jj-fnmin+1)=TimeAvDan(jc).DGAV(:,dgfind(1))./area;

            prof(j).ACu_W(:,jj)=getDGAVs('ACu_W',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); %getDGAVs finds correct data and divides by area
            
            
            prof(j).liq(:,jj)=getDGAVs('ALL_Q02',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); %all condensates
     end
    
    if cono==2332
            
%             col_str='ACu_W';
%             
%             area_str=strcat(colstr(1:3),'_A');
%             dgfind=findhead(area_str,dgstrDan(jc).dg);
%             area=TimeAvDan(jc).DGAV(:,dgfind(1));
%             areazeroes=find(area==0);
%             area(areazeroes)=1;
%             
%             dgfind=findhead(col_str,dgstrDan(jc).dg);
%             ACu_W(j).prof(:,jj-fnmin+1)=TimeAvDan(jc).DGAV(:,dgfind(1))./area;

            prof(j).ACu_W(:,jj)=getDGAVs('ACu_W',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); %getDGAVs finds correct data and divides by area
            
            
            prof(j).liq(:,jj)=getDGAVs('ALL_Q02',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); %all condensates
            prof(j).rain(:,jj)=getDGAVs('ALL_Q03',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            prof(j).snow(:,jj)=getDGAVs('ALL_Q04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            prof(j).graupel(:,jj)=getDGAVs('ALL_Q05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            prof(j).ice(:,jj)=getDGAVs('ALL_Q06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            
            Surf(j).MaxW(1:APARMS(2)+2,jj)=max(TwoDDan(jc).W,[],1);
            Surf(j).RainMR(1:APARMS(2)+2,jj)=TwoDDan(jc).Q(1,:,3);
            Surf(j).prec(1:APARMS(2),1:APARMS(1),jj)=surfdiag.prec;
            Surf(j).pmax(1:APARMS(2),1:APARMS(1),jj)=surfdiag.pmax;
            Surf(j).instant(1:APARMS(2),1:APARMS(1),jj)=surfdiag.instant;
            
            ztot=Radar(GridDan(jc),TwoDDan(jc));
        ztot=10*log10(ztot);
        for irad=1:size(TwoDDan(jc).W,2)-2
            ar=find(ztot(:,irad)>=10);
            if length(ar)>0
                echotops(irad,jj)=GridDan(jc).Z(ar(end));
            else
                echotops(irad,jj)=0;
            end
        end
        [minh ih]=min(abs(GridDan(jc).Z-3500));
        
        cappi(1:3,:,jj)=ztot(ih-1:ih+1,:);
        
            time(jj)=SerDan(jc).SER(end,1);
            
            
    end
    
    if strcmp(constr,'Wind Fields')==1
        
            times(jj)=SerDan(jc).SER(end,1); %store times of diags

            wind(j).W(:,:,jj)=TwoDDan(j).W;
            wind(j).V(:,:,jj)=TwoDDan(j).V;
            
%             dgfind=findhead('_A ',dgstrDan(jc).dg); %store areas of all the partitions
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('ALL_W ',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('ALL_V ',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('ALd_W ',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('ALu_W ',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
     end
     
     if strcmp(constr,'Vapour Fields')==1
        
            %times(jj)=SerDan(jc).SER(end,1); %store times of diags

            vap(j).v(:,:,jj)=TwoDDan(j).Q(:,:,1);
            pressure(j).p(:,:,jj)=TwoDDan(j).PP;
            potemp(j).p(:,:,jj)=TwoDDan(j).TH2; %store all potential temps
            
            icemr(j).i(:,:,jj)=1000*TwoDDan(jc).Q(:,:,6); %NOTE in g/kg
            icenc(j).i(:,:,jj)=TwoDDan(jc).Q(:,:,7);   
         
            
            
            times(jj)=SerDan(jc).SER(end,1); %store times of diags
            
%             dgfind=findhead('ALL_WQ',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 Fluxdiag(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end

%             dgfind=findhead('ALL_PSSUB',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 pssub(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end
%             
%             dgfind=findhead('ALL_PISUB',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 pisub(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end
%             
%             dgfind=findhead('ALL_PGSUB',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 pgsub(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end
%             
%             dgfind=findhead('ALL_PIMLT',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 pimlt(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end
%             
%             dgfind=findhead('ALL_PSMLT',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 psmlt(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end
%             
%             dgfind=findhead('ALL_PGMLT',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 pgmlt(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end
% %             
%             dgfind=findhead('ALL_WQ',dgstrDan(jc).dg); %store areas of all the partitions
%             for ifall=1:length(dgfind)
%                 Fluxdiag(j).dg(:,ifall,jj)=TimeAvDan(jc).DGAV(:,dgfind(ifall));
%             end
%             
            
            
%             dgfind=findhead('_A ',dgstrDan(jc).dg); %store areas of all the partitions
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('ALL_Q',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('ALu_Q',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('ALd_Q',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('CLu_Q',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgfind=findhead('CLd_Q',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
% 
%             dgfind=findhead('ACC_Q',dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
     end
     
     if strcmp(constr,'Temp')==1
        
            %times(jj)=SerDan(jc).SER(end,1); %store times of diags
            
            pressure(j).p(:,:,jj)=TwoDDan(j).PP;
            potemp(j).p(:,:,jj)=TwoDDan(j).TH2; %store all potential temps
            %gridsave(jj,j)=GridDan(j);
            
     end
     
    if strcmp(constr,'Ozone Fields')==1
        
            %times(jj)=SerDan(jc).SER(end,1); %store times of diags

            %potemp(j).p(:,:,jj)=TwoDDan(j).TH2; %store all potential temps
            %gridsave(jj,j)=GridDan(j);
            
            ozone(j).o(:,:,jj)=TwoDDan(jc).Q(:,:,14);
            
            pdat=1000*TwoDDan(jc).Q(:,:,2);
            for kk=3:6;
                pdat=pdat+1000*TwoDDan(jc).Q(:,:,kk); %in g/kg
            end;
            
            cond(j).c(:,:,jj)=pdat;
     end
     
     if strcmp(constr,'ice fields')==1
        
            icemr(j).i(:,:,jj)=1000*TwoDDan(jc).Q(:,:,6);
            icenc(j).i(:,:,jj)=TwoDDan(jc).Q(:,:,7);
            
            snowmr(j).i(:,:,jj)=1000*TwoDDan(jc).Q(:,:,4);
            snownc(j).i(:,:,jj)=TwoDDan(jc).Q(:,:,9);
            
            graupelmr(j).i(:,:,jj)=1000*TwoDDan(jc).Q(:,:,5);
            graupelnc(j).i(:,:,jj)=TwoDDan(jc).Q(:,:,8);
    
%         
     end
     
     if strcmp(constr,'snow graupel')==1
        
            %icemr(j).i(:,:,jj)=1000*TwoDDan(jc).Q(:,:,6);
            %icenc(j).i(:,:,jj)=TwoDDan(jc).Q(:,:,7);   
         
            snowmr(j).i(:,:,jj)=1000*TwoDDan(jc).Q(:,:,4);
            snownc(j).i(:,:,jj)=TwoDDan(jc).Q(:,:,9);
            
            graupelmr(j).i(:,:,jj)=1000*TwoDDan(jc).Q(:,:,5);
            graupelnc(j).i(:,:,jj)=TwoDDan(jc).Q(:,:,8);
    
%         
     end
     
     if strcmp(constr,'Sat MR & cond')==1
        
            %times(jj)=SerDan(jc).SER(end,1); %store times of diags

            %potemp(j).p(:,:,jj)=TwoDDan(j).TH2; %store all potential temps
            %gridsave(jj,j)=GridDan(j);
            
            switch j
            case 1
                pdat=satmr(2).s(:,:,jj);
                logcol(1)=1;
            case 2
                pdat=cond(2).c(:,:,jj);
                logcol(2)=0;
            end
            
            titlecb='';
            same=0;
            
            anim=1;
            
            %lims(1).l(1)=0.01; %log10(min(min(min(satmr(2).s))));
            %lims(2).l(1)=min(min(min(cond(2).c)));
            %lims(1).l(2)=log10(max(max(max(satmr(2).s))));
            %lims(2).l(2)=max(max(max(cond(2).c)));
                
     end
     
     if strcmp(constr,'Sat MR & cond')==1
        
            %times(jj)=SerDan(jc).SER(end,1); %store times of diags

            %potemp(j).p(:,:,jj)=TwoDDan(j).TH2; %store all potential temps
            %gridsave(jj,j)=GridDan(j);
            
            switch j
            case 1
                pdat=satmr(2).s(:,:,jj);
                logcol(1)=1;
            case 2
                pdat=cond(2).c(:,:,jj);
                logcol(2)=0;
            end
            
            titlecb='';
            same=0;
            
            anim=1;
            
            %lims(1).l(1)=0.01; %log10(min(min(min(satmr(2).s))));
            %lims(2).l(1)=min(min(min(cond(2).c)));
            %lims(1).l(2)=log10(max(max(max(satmr(2).s))));
            %lims(2).l(2)=max(max(max(cond(2).c)));
                
     end
     
     if strcmp(constr,'Vap MR & ice MR')==1
            
            time=mod(SER(end,1)/3600 + 19.67,24);
            timlab=num2str(round2(time,2));
            noplot=1;
            plotTimeHeightVap2;
            output=1;
            exdirA='c:/cygwin/home/user/runs/dmidamp_2/results/tot+iceNC/';
            
                
     end
     
     
     if strcmp(constr,'micro diags')==1
            
         noplot=1;
        output=0;
        
            %times(jj)=SerDan(jc).SER(end,1); %store times of diags

            %potemp(j).p(:,:,jj)=TwoDDan(j).TH2; %store all potential temps
            %gridsave(jj,j)=GridDan(j);
            
            switch j
            case 1
                pdat=satmr(2).s(:,:,jj);
                logcol(1)=1;
            case 2
                pdat=cond(2).c(:,:,jj);
                logcol(2)=0;
            end
            
            titlecb='';
            same=0;
            
            anim=1;
            
            %lims(1).l(1)=0.01; %log10(min(min(min(satmr(2).s))));
            %lims(2).l(1)=min(min(min(cond(2).c)));
            %lims(1).l(2)=log10(max(max(max(satmr(2).s))));
            %lims(2).l(2)=max(max(max(cond(2).c)));
                
     end
     
    
     
    if strcmp(constr,'Aerosol diags')==1
        %puts diags for all partitions of properties below into array diag for each time (puts into same positions as in DGAV)
        
        noplot=1;
        output=0;
        
        times(jj)=SerDan(jc).SER(end,1); %store times of diags
        
        dgfind=findhead('_A ',dgstrDan(jc).dg); %store areas of all the partitions
        diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
        
         dgfind=findhead('ALL_WQ',dgstrDan(jc).dg);
         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
         
         dgfind=findhead('ALu_WQ',dgstrDan(jc).dg);
         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
         
         dgfind=findhead('ALd_WQ',dgstrDan(jc).dg);
         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
         
         dgfind=findhead('CLu_WQ',dgstrDan(jc).dg);
         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
         dgfind=findhead('CLd_WQ',dgstrDan(jc).dg);
         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
         
         dgfind=findhead('ACC_WQ',dgstrDan(jc).dg);
         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
%         dgfind=findhead('aeroNP ',dgstrDan(jc).dg);
%         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
%         dgfind=findhead('NaeroF ',dgstrDan(jc).dg);
%         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
%         dgfind=findhead('Ndrops ',dgstrDan(jc).dg);
%         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
%          dgfind=findhead('QNloss',dgstrDan(jc).dg);
%          diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
% %         
% %         dgfind=findhead('Mloss ',dgstrDan(jc).dg);
% %         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
%         dgfind=findhead('maxS ',dgstrDan(jc).dg);
%         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
%         dgfind=findhead('Ndropmax ',dgstrDan(jc).dg);
%         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
% %         dgfind=findhead('MaxD ',dgstrDan(jc).dg);
% %         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
% %         
% %         dgfind=findhead('Ndropmax ',dgstrDan(jc).dg);
% %         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         
%         dgfind=findhead('Nact ',dgstrDan(jc).dg);
%         diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
        
%         for iq=15:54
%             dgstr=strcat('WQ',num2str(iq),' ');
%             dgfind=findhead(dgstr,dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%             
%             dgstr=strcat('Q',num2str(iq),' ');
%             dgfind=findhead(dgstr,dgstrDan(jc).dg);
%             diag(j).dg(:,dgfind,jj)=TimeAvDan(jc).DGAV(:,dgfind);
%         end
    end    
    
    if cono==23322
        ztot=Radar(GridDan(jc),TwoDDan(jc));
        ztot=10*log10(ztot);
        for irad=1:size(TwoDDan(jc).W,2)-2
            ar=find(ztot(:,irad)>=10);
            if length(ar)>0
                echotops(irad,jj)=GridDan(jc).Z(ar(end));
            else
                echotops(irad,jj)=0;
            end
        end
        [minh ih]=min(abs(GridDan(jc).Z-3500));
        
        cappi(1:3,:,jj)=ztot(ih-1:ih+1,:);
        time(jj)=SerDan(jc).SER(end,1);
    end
        
    
    if cono==20
        imrprof(j).prof(:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,6);
        incprof(j).prof(:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,7);
        gmrprof(j).prof(:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,5);
        gncprof(j).prof(:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,8);
        smrprof(j).prof(:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,4);
        sncprof(j).prof(:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,9);
        plotflag=1;
    
    end
    
    if cono==200
        maxLowTracer(j).prof(:,jj)=max(TwoDDan(jc).Q(:,:,10),[],2);
       noplot=1;
    
    end
    
    
    
    if cono==21
%         vlrprof(j).prof(1,:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,1); %vapour
%         vlrprof(j).prof(2,:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,2); %liquid
%         vlrprof(j).prof(3,:,jj-fnmin+1)=GridDan(jc).OLQBAR(:,3); %rain
   
        %plotflag=1;
        
            
        prof(j).vapmax(:,jj)=max(TwoDDan(jc).Q(:,:,1),[],2);     %max vapour
        prof(j).vapmin(:,jj)=min(TwoDDan(jc).Q(:,:,1),[],2);     %min vapour
        prof(j).icemax(:,jj)=max(TwoDDan(jc).Q(:,:,6),[],2);     %max ice MR 
        prof(j).snowmax(:,jj)=max(TwoDDan(jc).Q(:,:,4),[],2);     %max snow MR
        prof(j).graupmax(:,jj)=max(TwoDDan(jc).Q(:,:,5),[],2);     %max graupel MR
        prof(j).iceNCmax(:,jj)=max(TwoDDan(jc).Q(:,:,7),[],2);     %max ice NC 
        prof(j).snowNCmax(:,jj)=max(TwoDDan(jc).Q(:,:,9),[],2);     %max snow NC
        prof(j).graupNCmax(:,jj)=max(TwoDDan(jc).Q(:,:,8),[],2);     %max graupel NC
        prof(j).midTRmax(:,jj)=max(TwoDDan(jc).Q(:,:,11),[],2);     %max mid tracer
        prof(j).uppTRmax(:,jj)=max(TwoDDan(jc).Q(:,:,12),[],2);     %max upper tracer

        prof(j).vapmean(:,jj)=mean(TwoDDan(jc).Q(:,:,1),2);     %mean vapour MR    
        prof(j).icemean(:,jj)=mean(TwoDDan(jc).Q(:,:,6),2);     %mean ice MR
        prof(j).snowmean(:,jj)=mean(TwoDDan(jc).Q(:,:,4),2);     %mean snow MR
        prof(j).graupmean(:,jj)=mean(TwoDDan(jc).Q(:,:,5),2);     %mean graupel MR
        prof(j).iceNCmean(:,jj)=mean(TwoDDan(jc).Q(:,:,7),2);     %mean ice NC
        prof(j).snowNCmean(:,jj)=mean(TwoDDan(jc).Q(:,:,9),2);     %mean snow NC
        prof(j).graupNCmean(:,jj)=mean(TwoDDan(jc).Q(:,:,8),2);     %mean graupel NC
        prof(j).midTRmean(:,jj)=mean(TwoDDan(jc).Q(:,:,11),2);     %mean ice MR
        prof(j).uppTRmean(:,jj)=mean(TwoDDan(jc).Q(:,:,12),2);     %mean ice MR

        noplot=1;
        
    end
    
    if cono==22
        ALLDTH(j).prof(1,:,jj-fnmin+1)=TimeAvDan(jc).DGAV; %vapour
       
        plotflag=2;
    
    end
    
    if strcmp(lower(constr),lower('Ice Proc diags ALL'))==1
        %icediags_5thSept_2005_2;
        icediags_5thSept_2005_32;
    end
    
    
    if strcmp(constr,'Ice Proc diags')==1
        clear dgs;
%         dgs{1}='PIMLT';
% 		dgs{2}='PSAUT';
% 		dgs{3}='PSACI';
% 		dgs{4}='PRACI_S';
% 		dgs{5}='PGACI';
% 		dgs{6}='PRACI_G';
% 		dgs{7}='PIHAL';
% 		dgs{8}='PIPRM';
% 		dgs{9}='PICNT';
% 		dgs{10}='PIDEP';
% 		dgs{11}='PIACW';
% 		dgs{12}='PIFRW';
%         dgs{13}='RSAUT';
%         dgs{14}='RIACI';
        
        dgs{1}='PISUB';
		dgs{2}='PSDEP';
		dgs{3}='PIACR_S';
		dgs{4}='PSACW';
		dgs{5}='PSSUB';
		dgs{6}='PGACS';
		dgs{7}='PRACS';
		dgs{8}='PGAUT';
		dgs{9}='PSMLT';
		dgs{10}='RSBRK';
		dgs{11}='RIACR_S';
		dgs{12}='RGACS';
        dgs{13}='RSACR';
        dgs{14}='RSACS';
        
        for iprc=1:length(dgs);
            nam=strcat('ALL_',dgs{iprc});
            icediag(j).i(:,jj,iprc)=getDGAVs(nam,dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
        end
        
    end
    
    if strcmp(constr,'Ice NC anim')==1
        switch j
        case 1
            %pdat=icenc(1).i(:,:,jj);
            pdat=TwoDDan(1).Q(:,:,7);
            minZovr(j)=14e3;
            maxZovr(j)=23e3;
            titles(1).tit='Ice Number Concentration (#/kg)';
            titlecb='';
        case 2
            pdat=sum(TwoDDan(1).Q(:,:,2:5),3)*1000;
            minZovr(j)=0e3;
            maxZovr(j)=30e3;
            titles(2).tit='Total Condensate (g/kg)';
            TwoDDan(2)=TwoDDan(1);
            GridDan(2)=GridDan(1);
            titlecb='';
        end
    end
    
    if strcmp(constr,'Ice Proc diags2')==1
        
            %ice MR
            icediag2(j).i(:,jj,1)=getDGAVs('ALu_DQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            icediag2(j).i(:,jj,2)=getDGAVs('ALu_WQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,3)=getDGAVs('ALu_Q06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            
            %ice NC
            icediag2(j).i(:,jj,4)=getDGAVs('ALu_DQ07',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,5)=getDGAVs('ALu_WQ07',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            icediag2(j).i(:,jj,6)=getDGAVs('ALu_Q07',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            
            %snow MR
            icediag2(j).i(:,jj,7)=getDGAVs('ALu_DQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            icediag2(j).i(:,jj,8)=getDGAVs('ALu_WQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,9)=getDGAVs('ALu_Q04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            
            %snow NC
            icediag2(j).i(:,jj,10)=getDGAVs('ALu_DQ09',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,11)=getDGAVs('ALu_WQ09',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            icediag2(j).i(:,jj,12)=getDGAVs('ALu_Q09',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            
            %means from TwoD
            %icediag2(j).i(:,jj,13)=mean(TwoD.Q(:,:,6),2);
            %icediag2(j).i(:,jj,14)=mean(TwoD.Q(:,:,7),2);
            %icediag2(j).i(:,jj,15)=mean(TwoD.Q(:,:,4),2);
            %icediag2(j).i(:,jj,16)=mean(TwoD.Q(:,:,9),2);
 
        
    end
    
    if strcmp(constr,'Ice Proc diags3')==1
        
            %ice MR
            icediag2(j).i(:,jj,1)=getDGAVs('ALL_WQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            icediag2(j).i(:,jj,2)=getDGAVs('ALL_Q06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,3)=getDGAVs('ALL_FQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            
            %snow MR
            icediag2(j).i(:,jj,4)=getDGAVs('ALL_WQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            icediag2(j).i(:,jj,5)=getDGAVs('ALL_Q04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,6)=getDGAVs('ALL_FQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            
            %graupel MR
            icediag2(j).i(:,jj,7)=getDGAVs('ALL_WQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,8)=getDGAVs('ALL_Q05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            icediag2(j).i(:,jj,9)=getDGAVs('ALL_FQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            
            %vapour
            icediag2(j).i(:,jj,10)=getDGAVs('ALL_WQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
            icediag2(j).i(:,jj,11)=getDGAVs('ALL_Q01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            %icediag2(j).i(:,jj,12)=getDGAVs('ALu_FQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            
            %means from TwoD
            %icediag2(j).i(:,jj,13)=mean(TwoD.Q(:,:,6),2);
            %icediag2(j).i(:,jj,14)=mean(TwoD.Q(:,:,7),2);
            %icediag2(j).i(:,jj,15)=mean(TwoD.Q(:,:,4),2);
            %icediag2(j).i(:,jj,16)=mean(TwoD.Q(:,:,9),2);
                
            time(j,jj)=SerDan(j).SER(end,1);
        
    end
    
    if strcmp(constr,'Ice prctiles')==1
        id=1;
        prcs=[0:5:100];
        vap_prctiles(j).v(:,jj,id:id+length(prcs)-1)=(prctile(TwoDDan(j).Q(:,:,1)',prcs))';
        time(j,jj)=SerDan(j).SER(end,1);
    end   
        
        
    
    if strcmp(constr,'Ice Proc diags4')==1
        
            id=1;
         
            
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQSG01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQSG01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_FQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_FQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             
%             
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQSG04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQSG04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_FQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_FQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQSG05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQSG05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_FQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_FQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_WQSG06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin); 
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_WQSG06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_FQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_FQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_A',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
%             icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_A',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
%             id=id+1;
            
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_DQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_DQ01',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_DQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_DQ04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_DQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_DQ05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALu_DQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            icediag4(j).i(:,jj,id)=getDGAVsNoArea('ALd_DQ06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,jj,fnmin);
            id=id+1;
            
%             
            %NOTE getDGAVs divides by the areas!!!
            %getDGAVsNoArea doesn't
            
            %means from TwoD
            %icediag2(j).i(:,jj,13)=mean(TwoD.Q(:,:,6),2);
            %icediag2(j).i(:,jj,14)=mean(TwoD.Q(:,:,7),2);
            %icediag2(j).i(:,jj,15)=mean(TwoD.Q(:,:,4),2);
            %icediag2(j).i(:,jj,16)=mean(TwoD.Q(:,:,9),2);
 
        
    end
    
    
    
    
    if cono==113
        pdat=TwoDDan(jc).W;
        titlecb='W (m/s)';
    end;
    
    if cono==114
        pdat=TwoDDan(jc).TH1;
        pdat(1,:)=0;
        titlecb='Theta Pert. (K)';
    end;
   
    
    
    if cono==16;
        titlecb='Total RSGI MR (/kg)';
        pdat=TwoDDan(jc).Q(:,:,3);
        for kk=4:6;
            pdat=pdat+TwoDDan(jc).Q(:,:,kk);
        end;
        
    end;
    if cono==17;
        titlecb='Total NC (/kg)';
        pdat=TwoDDan(jc).Q(:,:,6);
        for kk=14:15;
            pdat=pdat+TwoDDan(jc).Q(:,:,kk-7);
        end;
    end;
    
    if cono==18;
            titlecb='Tracer MR (kg/kg)';
            pdat=TwoDDan(jc).Q(:,:,10);
    end;
    
    if cono==19;
            titlecb='Radar reflectivity dBz';
            TwoD=TwoDDan(jc);
            ztot=radar(GridDan(jc),TwoDDan(jc));
            pdat=10*log10(ztot);
            maxRad(j).prof(:,jj)=max(pdat,[],2); %max reflectivity profile
    end;
    
    if length(logcol)>=j
        logcolj=logcol(j);
    else
        logcolj=logcol(1);
    end
    
    if logcolj==1 & noplot==0
        pdat=log10(pdat);
    end
    
   
    if noplot==0
        if maxZovr(j)>-1
                iz=findheight(GridDan(jc).Z,maxZovr(j));
        else
                iz=size(GridDan(jc).Z,1);
        end
        if minZovr(j)>-1
                izmin=findheight(GridDan(jc).Z,minZovr(j));
        else
                izmin=1;
        end
        maxims(j)=max(max(pdat(izmin:iz,:)));
        minims(j)=min(min(pdat(izmin:iz,:)));
        if useozmin==1 & cono==266 & qf==133; minims(j)=mq*x0 + kq; end
    end
        
    
    %pdat(1,1)=maxval;
    if noplot==0
    
    if plotflag==0
        if contours==1
           % contourf(XX./1000,ZZ./1000,pdat);
           [hcon ccon]=contour(XX./1000,ZZ./1000,pdat,[300:5:390 400:10:500]);clabel(hcon,ccon)
        else
		    pcolor(XX./1000,ZZ./1000,pdat);
            shading interp;
        end
        
        maxval(jj)=max(maxims);
	    minval(jj)=min(minims);
	
    if samesize==1
        maxXX=max(xxmax);
	    maxZZ=max(zzmax);
    else
        maxXX=max(xxmax(jc));
	    maxZZ=max(zzmax(jc));
    end


	if maxZovr(j)>-1
	    maxZZ=maxZovr(j); %override maxZZ
	end
    
    if minZovr(j)>-1
	    minZZ=minZovr(j); %override maxZZ
    else
        minZZ=0;
    end
    
    if limcol==1
         minval(jj)=limcols(1);
         maxval(jj)=limcols(2);
    end

    
        
        if same~=1
            izmin=findheight(GridDan(jc).Z,minZZ);
            izmax=findheight(GridDan(jc).Z,maxZZ);
            ixmin=findheight(GridDan(jc).Y1,-maxXX);
            ixmax=findheight(GridDan(jc).Y1,maxXX);    
            
            lims(j).l(1)=minALL(pdat(izmin:izmax,ixmin:ixmax));
            lims(j).l(2)=maxALL(pdat(izmin:izmax,ixmin:ixmax));
        end
        
            %axis equal;
		if (cono==16|cono==19|cono==5|cono==26) %& contours==0
            d=max(max(TwoDDan(jc).Q(:,:,3)));
            hold on;
            contour(XX./1000,ZZ./1000,TwoDDan(jc).Q(:,:,3),[d/20 d/20],'w');
            contour(XX./1000,ZZ./1000,TwoDDan(jc).Q(:,:,3),[d/1.1 d/1.1],'w');
		end
		
		
		if streamflag==1
            VW=TwoDDan(jc).V;
            WW=TwoDDan(jc).W;
            size(VW);
            %hss=streamslice(XX./1000,ZZ./1000,VW,WW,ans(1)/150);
            hss=streamslice(XX./1000,ZZ./1000,VW,WW,1);
            set(hss,'Color','black');
		end
        
        if vectorf==1
             if maxZovr>1
                iz=find(Grid.Z>maxZovr(j));
                iz=iz(1);
            else
                iz=size(Grid.Z,1);
            end
            
            sqy=size(Grid.Y1,1)/25;
            sqz=round(iz/15);
            
            sf=max(max(TwoDDan(jc).V))/max(max(TwoDDan(jc).W));
           
%             clear meanV;
%             for ih=1:size(Grid.Z,1)
%                 meanV(ih,:)=TwoDDan(jc).V(ih,:)-GridDan(jc).VBAR(ih);
%             end
            hold on;
            %for iq=1:sqz:iz
             %   qdat(iq,:)=TwoDDan(jc).V(iq,:)-Grid.VBAR(iq);
             %end
            quiverDan(GridDan(jc).Y1(1:sqy:end)./1000,GridDan(jc).Z(1:sqz:end)./1000,TwoDDan(jc).V(1:sqz:end,1:sqy:end),TwoDDan(jc).W(1:sqz:end,1:sqy:end),'w');
        end
		
		%title(tit(j).tit,'fontsize',7);
		%hc((jj-1)*nplots+j).h=colorbar;
		
		%if j==5|j==11;
		xlabel(texlabel('Horizontal Displacement(km)'),'fontsize',fsize);
		ylabel( texlabel('Height(km)'),'fontsize',fsize);
		
		%end;
       
		
        if ititle==1
            title(titles(j).tit,'fontsize',fsize);
        end
    
	else % if plotflag==0
        xdat=Grid.Y1(1:sb:sb*nb)./1000;
        hold on;
        %plot(xdat,bin);
        %plot(xdat,bin,'kx');
	end %if plotflag==0

end %noplot==0

 end; %subplot end for3

 if noplot==0
    axes('Position',[0 0 1 1],'Visible','off');
    strcat('TIME = ',num2str( SerDan(idir).SER(end,1)/3600 ),' hours');
    if itext==1
        text(.4,.985,strcat('TIME = ',num2str( SerDan(idir).SER(end,1)/3600 ),' hours'));
    end
 end

end;  %figure end for2




if (plotflag==0 & noplot==0)
   'here2'     
	

	for j=1:nplots
        
        maxval(jj)=max(maxims);
	minval(jj)=min(minims);
	
    if samesize==1
        maxXX=max(xxmax);
	    maxZZ=max(zzmax);
    else
        maxXX=max(xxmax(jc));
	    maxZZ=max(zzmax(jc));
    end


	if maxZovr(j)>-1
	    maxZZ=maxZovr(j); %override maxZZ
	end
    
    if minZovr(j)>-1
	    minZZ=minZovr(j); %override maxZZ
    else
        minZZ=0;
    end
    
    if limcol==1
         minval(jj)=limcols(1);
         maxval(jj)=limcols(2);
    end
    
    
        
        if zlim(j)==1
            axis(hs((jj-1)*nplots+j).h,[-maxXX/1000 maxXX/1000 minZZ/1000 maxZZ/1000]);
        end
        set(hs((jj-1)*nplots+j).h,'fontsize',fsize);
        %caxis(hs((jj-1)*nplots+j).h,[0 maxval(jj-1)]);
        
        
        
        
        
        if contours==0
            if cono==19
                caxis(hs((jj-1)*nplots+j).h,[0 60]);   
            elseif same==1
                if minval(jj)>=maxval(jj); minval(jj)=maxval(jj)-1;end
                caxis(hs((jj-1)*nplots+j).h,[minval(jj) maxval(jj)]);
            elseif length(lims)>1
                if lims(j).l(2)>lims(j).l(1)
                    caxis(hs((jj-1)*nplots+j).h,[lims(j).l(1) lims(j).l(2)]);
                end
            end
            
            hc=colorbar('peer',hs((jj-1)*nplots+j).h);
            
            if length(logcol)>=j
                logcolj=logcol(j);
            else
                logcolj=logcol(1);
            end
    
            if logcolj==1
                clear ctickstr
                ctick=get(hc,'ytick');
                for j=1:length(ctick)
                    %te=strcat('10^','{',num2str(ctick(j)),'}' );
                    te=num2str(10^ctick(j),3);
                    ctickstr(j,1:length(te))=te;
                end
                
                set(hc,'yticklabel','');
                add=ctick(end)/50;
                set(hff,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
                text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );
			end
            
            %set(hc,'fontsize',fsize-6);
            
            
            set(hc,'fontsize',fsize);
            
        end
        
        if cono==266 & qf==133
            clear cbtick;
            cbvals=get(hc,'ytick');
            cbdiff=cbvals(2)-cbvals(1);
            cbvals=[cbvals(1)-cbdiff:cbdiff/2:cbvals(end)+cbdiff];
            acb=find(cbvals<y1);
            tvals=(cbvals(acb)-kq)/mq;
            acb=find(cbvals>y1);
            tvals(end+1:end+length(acb))=exp(cbvals(acb));
            for icb=1:length(tvals)
                te=num2str(tvals(icb),2);
                cbtick(icb,1:length(te))=te;
            end
            set(hc,'ytick',cbvals);
            set(hc,'yticklabels',cbtick);
        end
            
            
	end;
   
    
    
    text(maxXX*1.1/1000,-maxZZ*0.1/1000,titlecb,'fontsize',fsize);
    
     if anim==1
		b(jj-fnmin+1)=getframe;
        close(gcf);
     end
    
        
        
        
end %(plotflag==0 & noplot==0)

        if output==1
			for k=1:kmax
                textdataDan(1).text=strcat('TIME=',num2str(SER(end,1)/3600));
				if field==0
                    exname=strcat(exdirA,int2str(jj),'-',int2str(k),'-',textdataDan(1).text,'-',int2str(cono),'-',constr,'.emf'); 
				else
                    exname=strcat(exdirA,direcDan(k).dir(12:end),int2str(k),'-',textdataDan(1).text,'-','.jpg');
				end
				
				if exist('hff'); gcf=hff;  end %hh((jj-1)*kmax+k).h;
				set(gcf,'paperpositionmode','auto');
				%print(gcf,'-djpeg','-r350',exname);
                print(gcf,'-dmeta',exname);
                
                %print(gcf,'-dmeta',exname);
                
				close(gcf);
			end
        end %if output==1      
        %text(maxXX*0.9,-maxZZ/10,titlecb,'units','centimeters','fontsize',fsize);
        
    
end; %fnmax (jj) end for1

if (plotflag==1 & noplot==0)
    xdat=Grid.Y1(1:sb:sb*nb)./1000;
    figure;
    for jj=1:nplots
        hs(jj).h=subplot(a,b,jj); %subplotting
        plot(xdat,bin(jj).b);
        hold on;
        plot(xdat,bin(jj).b,'kx')
        title(direcDan(jc).dir);
	end
    
    figure;
    for jj=1:nplots
        hs(jj).h=subplot(a,b,jj); %subplotting
        plot(time,bintime(jj).b);
        hold on;
        plot(time,bintime(jj).b,'kx');
        title(direcDan(jc).dir);
	end
end


break;
for jj=2:fnmax
    for j=1:nplots
        %caxis(hs((jj-1)*nplots+j).h,[0 max(maxval)]);
        caxis(hs((jj-1)*nplots+j).h,[0 maxval(jj)]);
        colorbar('peer',hs((jj-1)*nplots+j).h);
        axes('Position',[0 0 1 1],'Visible','off');
        text(.4,.99,textdataDan(1).text);
    end
   % bb(jj-1)=getframe(hh(jj));
end


