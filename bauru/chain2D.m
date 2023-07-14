%imports all diagnostic files up to a max and plots 2-d plots of each
titsbauru;
clear b hh;
fnmin=input('Enter the first diagnostic file number, with directories in direcDan(i).dir: ');
fnmax=input('Enter the maximum diagnostic file number, with directories in direcDan(i).dir: ');
con=input('Enter col no.: ');
cono=con;
nplots=input('Enter the no. of directories to import :');

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.1];


for jj=fnmin:fnmax;  %for1
    if jj<10;
        fn=strcat('0',int2str(jj));
    else;
        fn=int2str(jj);
    end;

    for i=1:nplots;
        FileName=strcat('c:\cygwin\home\user\',direcDan(i).dir,'\RUN0001.DG00',fn);
        DIAG2_3_v2;
        GridDan(i)=Grid;
        TwoDDan(i)=TwoD;
        SerDan(i).SER=SER;
        textdataDan(i).text=textdat;
    end;


clear maxims;   

kmax=fix((nplots-0.1)/6)+1;
for k=1:kmax     %for2
    hh((jj-2)*kmax+k).h=figure('name',strcat('2d plots with same scales of column:',int2str(cono)),'Position',posit);
    
    jmax=6;
    if nplots-k*6<0 
        jmax=rem(nplots,6)+(k-1)*6;
    end;
    
 for j=(k-1)*6+1:jmax*k;  %for3
    plotnum=rem(j,6);
    if plotnum==0;
        plotnum=6;
    end
    
    hs((jj-2)*nplots+j).h=subplot(3,2,plotnum);
    axis equal;
    [ZZ XX]=meshgrid(GridDan(j).Z,GridDan(j).Y1);
    ZZ=ZZ';
    XX=XX';
    xxmax(j)=max(GridDan(j).Y1)/1000;
    zzmax(j)=max(GridDan(j).Z)/1000;
    
    %for j=1:nplots
    
    if cono==3
        pdat=TwoDDan(j).U;
        titlecb='U (m/s)';
    end;
    if cono==4
        pdat=TwoDDan(j).V;
        titlecb='V (m/s)';
    end;
    if cono==5
        pdat=TwoDDan(j).W;
        titlecb='W (m/s)';
    end;
    if cono==6
        pdat=TwoDDan(j).TH;
        titlecb='Theta (K)';
    end;
   
    if cono>=7 & cono<=15
        if cono==9;
            titlecb='Rain Water MR (kg/kg)';
        end;
        if cono==10;
            titlecb='Snow MR (kg/kg)';
        end;
        if cono==11;
            titlecb='Graupel MR (kg/kg)';
        end;
        if cono==12;
            titlecb='Ice MR (kg/kg)';
        end;
        if cono==13;
            titlecb='Ice NC (/kg)';
        end;
        if cono==14;
            titlecb='Graupel NC (/kg)';
        end;
        if cono==15;
            titlecb='Snow NC (/kg)';
        end;
        
    pdat=TwoDDan(j).Q(:,:,cono-6);
    
        if cono==7;
            titlecb='Vapour MR (kg/kg)';
            pdat(:,1)=0;
        end;
        if cono==8;
            titlecb='Liquid MR (kg/kg)';
            pdat(:,1)=0;
        end;
    end;
    
     if cono==16;
        titlecb='Total of RSGI MRs (kg/kg)';
        pdat=TwoDDan(j).Q(:,:,3);
        for kk=11:13;
            pdat=pdat+TwoDDan(j).Q(:,:,kk-7);
        end;
    end;
    if cono==17;
        titlecb='Total NC (/kg)';
        pdat=TwoDDan(j).Q(:,:,6);
        for kk=14:15;
            pdat=pdat+TwoDDan(j).Q(:,:,kk-7);
        end;
    end;
    
    if cono==18;
            titlecb='Tracer MR (kg/kg)';
            pdat=TwoDDan(j).Q(:,:,9);
    end;
    
    if cono==19;
            titlecb='Radar reflectivity dBz';
            TwoD=TwoDDan(j);
            radar;
            pdat=10*log10(ztot);
    end;
        
 
    maxims(j)=max(max(pdat));
        
    
    %pdat(1,1)=maxval;
    pcolor(XX./1000,ZZ./1000,pdat);
    if cono==16|cono==19
        d=max(max(TwoDDan(j).Q(:,:,3)));
        hold on;
        contour(XX./1000,ZZ./1000,TwoDDan(j).Q(:,:,3),[d/20 d/20],'k');
    end
    
    VW=TwoDDan(j).V;
    WW=TwoDDan(j).W;
    size(VW);
    %hss=streamslice(XX,ZZ,VW,WW,ans(1)/300);
    %set(hss,'Color','black');
    
    title(tit(j).tit,'fontsize',7);
    %hc((jj-1)*nplots+j).h=colorbar;
    shading interp;
    if j==5;
    xlabel(texlabel('Horizontal Displacement(km)'),'fontsize',7);
    ylabel( texlabel('Height(km)'),'fontsize',7);
    
    end;
    text(5.5,-0.3,titlecb,'units','centimeters','fontsize',7);
 %   b(jj-1)=getframe;
 
 end; %subplot end for3
 
axes('Position',[0 0 1 1],'Visible','off');
        text(.4,.98,'Time = 67 mins');
end;  %figure end for2
maxval(jj-1)=max(maxims);
maxXX=max(xxmax);
maxZZ=max(zzmax);

for j=1:nplots
    axis(hs((jj-2)*nplots+j).h,[-maxXX maxXX 0 maxZZ]);
    set(hs((jj-2)*nplots+j).h,'fontsize',7);
    caxis(hs((jj-2)*nplots+j).h,[0 maxval(jj-1)]);
    hc=colorbar('peer',hs((jj-2)*nplots+j).h);
    set(hc,'fontsize',7);
        
end;

for k=1:kmax
exname=strcat('c:\matlab6p1\work\',exdir,textdataDan(1).text,'-',int2str(k));    
gcf=hh((jj-2)*kmax+k).h
print(gcf,'-djpeg','-r150',exname);
close(gcf);
end

end; %fnmax (jj) end for1

break;
for jj=2:fnmax
    for j=1:nplots
        %caxis(hs((jj-1)*nplots+j).h,[0 max(maxval)]);
        caxis(hs((jj-2)*nplots+j).h,[0 maxval(jj-1)]);
        colorbar('peer',hs((jj-2)*nplots+j).h);
        axes('Position',[0 0 1 1],'Visible','off');
        text(.4,.97,textdataDan(1).text);
    end
   % bb(jj-1)=getframe(hh(jj));
end


