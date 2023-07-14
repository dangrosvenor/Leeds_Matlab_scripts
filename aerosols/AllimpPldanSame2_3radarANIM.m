%imports several diagnostic files, plots 2d, saves as jpeg and closes down figure


sorthead=0; %flag to switch off sorting of headers to save time
same=1; %flag for giving all plots the same colour scale (=1)
field=0; %flag for baurufield to add directories to filename
fsize=15; %fontsize
plotflag=0;
maxZovr=20000; %set to zero for top of domain
minZovr=500;
noplot=0; %flag to stop figures appearing (=1)
output=1; %flag for output to file (=1)
contoura=[0 0]; %array to tell which plot to use contours on
samesize=0; %flag to plot all 2d plots on the same size domain for comparing runs with different domain sizes
anim=1; %flag to make animation frames in animframe(jj-fmin+1) 

limcol=0;
limcols=[7350 7650];

itext=0;
ititle=0;

logcol=0; %flag to plot log of colour values

useozmin=1; %flag to use value below as min value on colour scale for ozone plots
ozmin=5.3e-8;

nb=20;


iread=1; %flag to say whether to read in files or just plot the one that's already in memory
streamflag=0; %flag for putting streamlines of wind on
vectorf=1;

%exdir='field\etasound\30.01-16\etavars\Radar\'; %backslash at end
exdir='c:\matlabR12\work\bauru/tracersjan2005/aerosols/Nenes_HeightVar/droplets'; % \ at end !
exdir='c:/cygwin/home/login/runs/aeroruns/aeroqadv2_4_deplete/results/droplets/'
exdir='c:\matlabR12\work\bauru/tracersjan2005/egugraphs/';
exdir='c:\upload/egugraphs/anim/force+0_3th3qv';

mov = avifile('c:\upload/egugraphs/anim/force+0_3th3qv/aerosolDep_1_FPS.avi','fps',1);

if strcmp(exdir(end),'/')==0; exdir(end+1)='/'; end

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


if cono==266
    qf=input('Enter Q field no. : ');
    qfstr=int2str(qf);
else
    qfstr='';
end

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];
%titsbaurufull;

sortheaderdgs;
%ndth=strmatch('ALL_DTH',dgstr); %column no. in TimeAv.DGAV for latent heat profiles due to microphysics




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
        
        titles(i).tit=direcDan(idir).dir;
       
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
        hh((jj-1)*kmax+k).h=figure('name',strcat('2d plots with same scales of column:',int2str(cono),' :',qfstr));
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
        hs((jj-2)*nplots+j).h=subplot(b,a,plotnum); %subplotting
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
            titlecb='Ice NC (/kg)';
        end;
        if cono==8;
            titlecb='Graupel NC (/kg)';
        end;
        if cono==9;
            titlecb='Snow NC (/kg)';
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
    
    if cono==26;
        titlecb='Total LRSGI MR (/kg)';
        pdat=TwoDDan(jc).Q(:,:,2);
        for kk=3:6;
            pdat=pdat+TwoDDan(jc).Q(:,:,kk);
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
    
    
    if logcol==1
        pdat=log(pdat);
    end
    
   
    if noplot==0
        if maxZovr>1
                iz=find(GridDan(jc).Z>maxZovr);
                iz=iz(1);
        else
                iz=size(GridDan(jc).Z,1);
        end
        if minZovr>1
                izmin=find(GridDan(jc).Z<minZovr);
                izmin=izmin(end);
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
                iz=find(Grid.Z>maxZovr);
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
            %quiverDan(GridDan(jc).Y1(1:sqy:end)./1000,GridDan(jc).Z(1:sqz:iz)./1000,TwoDDan(jc).V(1:sqz:iz,1:sqy:end),TwoDDan(jc).W(1:sqz:iz,1:sqy:end),'w');
        end
		
		%title(tit(j).tit,'fontsize',7);
		%hc((jj-1)*nplots+j).h=colorbar;
		
		%if j==5|j==11;
		xlabel(texlabel('Horizontal Displacement(km)'),'fontsize',fsize);
		ylabel( texlabel('Height(km)'),'fontsize',fsize);
		
		%end;
		%   b(jj-1)=getframe;
		
        if ititle==1
            title(titles(j).tit);
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
	maxval(jj-1)=max(maxims);
	minval(jj-1)=min(minims);
	
    if samesize==1
        maxXX=max(xxmax);
	    maxZZ=max(zzmax);
    else
        maxXX=max(xxmax(jc));
	    maxZZ=max(zzmax(jc));
    end


	if maxZovr>1
	    maxZZ=maxZovr; %override maxZZ
	end
    
    if minZovr>1
	    minZZ=minZovr; %override maxZZ
	end
    
    if limcol==1
         minval(jj-1)=limcols(1);
         maxval(jj-1)=limcols(2);
    end

	for j=1:nplots
        
        axis(hs((jj-2)*nplots+j).h,[-maxXX/1000 maxXX/1000 minZZ/1000 maxZZ/1000]);
        set(hs((jj-2)*nplots+j).h,'fontsize',fsize);
        %caxis(hs((jj-2)*nplots+j).h,[0 maxval(jj-1)]);
        
        if contours==0
            if cono==19
                caxis(hs((jj-2)*nplots+j).h,[0 65]);   
            elseif same==1
                if minval(jj-1)>=maxval(jj-1); minval(jj-1)=maxval(jj-1)-1;end
                caxis(hs((jj-2)*nplots+j).h,[minval(jj-1) maxval(jj-1)]);
            end
            
            hc=colorbar('peer',hs((jj-2)*nplots+j).h);
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
    
    if anim==1
        animfr(jj)=getframe(hh(jj).h);
        mov = addframe(mov,animfr(jj)); 
    end
   
    
    
    
        if output==1
			for k=1:kmax
                textdataDan(1).text=strcat('TIME=',num2str(SER(end,1)/3600));
				if field==0
                    exname=strcat(exdir,textdataDan(1).text,'-',int2str(k),'-',int2str(jj),'-',int2str(cono),'-',constr,'.jpg'); 
				else
                    exname=strcat(exdir,direcDan(k).dir(12:end),textdataDan(1).text,'-',int2str(k),'.jpg');
				end
				
				gcf=hh((jj-1)*kmax+k).h;
				set(gcf,'paperpositionmode','auto');
				print(gcf,'-djpeg','-r350',exname);
                
                %print(gcf,'-dmeta',exname);
                
				close(gcf);
			end
        end %if output==1      
        %text(maxXX*0.9,-maxZZ/10,titlecb,'units','centimeters','fontsize',fsize);
        text(maxXX*1.1/1000,-maxZZ*0.1/1000,titlecb,'fontsize',fsize);
        
end %(plotflag==0 & noplot==0)
    
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

mov = close(mov);

break;
for jj=2:fnmax
    for j=1:nplots
        %caxis(hs((jj-1)*nplots+j).h,[0 max(maxval)]);
        caxis(hs((jj-2)*nplots+j).h,[0 maxval(jj-1)]);
        colorbar('peer',hs((jj-2)*nplots+j).h);
        axes('Position',[0 0 1 1],'Visible','off');
        text(.4,.99,textdataDan(1).text);
    end
   % bb(jj-1)=getframe(hh(jj));
end


