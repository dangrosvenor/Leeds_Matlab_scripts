%imports all diagnostic files up to a max and plots 2-d plots of each

clear b;
fnmax=input('Enter the maximum diagnostic file number, with directories in direcDan(i).dir: ');
con=input('Enter col no.: ');
cono=con;
nplots=input('Enter the no. of directories to import :');

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];


for jj=2:fnmax;
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
    
for k=1:fix((nplots-0.1)/6)+1
    hh(jj)=figure('name',strcat('2d plots with same scales of column:',int2str(cono)),'Position',posit);
    
    jmax=6;
    if nplots-k*6<0 
        jmax=rem(nplots,6)+(k-1)*6;
    end;
    
for j=(k-1)*6+1:jmax;
    plotnum=rem(j,6);
    if plotnum==0;
        plotnum=6;
    end
    
    hs((jj-1)*nplots+j).h=subplot(1,1,plotnum);
    axis equal;
    [ZZ XX]=meshgrid(GridDan(j).Z,GridDan(j).Y1);
    ZZ=ZZ';
    XX=XX';
    xxmax(j)=max(GridDan(j).Y1);
    zzmax(j)=max(GridDan(j).Z);
    
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
        titlecb='Total MR (/kg)';
        pdat=TwoDDan(j).Q(:,:,2);
        for kk=10:12;
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
        
 
    maxims(j)=max(max(pdat));
        
    
    %pdat(1,1)=maxval;
    pcolor(XX,ZZ,pdat);
    
    VW=TwoDDan(j).V;
    WW=TwoDDan(j).W;
    size(VW);
    hss=streamslice(XX,ZZ,VW,WW,ans(1)/300);
    set(hss,'Color','black');
    
    title(direcDan(j).dir);
    colorbar;
    shading interp;
    if j==1;
    xlabel(texlabel('Horizontal Displacement(m)'));
    ylabel( texlabel('Height(m)') );
    text(10,-0.5,titlecb,'units','centimeters');
    end;
 %   b(jj-1)=getframe;
 
end; %subplot end
maxval(jj)=max(maxims);
maxXX=max(xxmax);
maxZZ=max(zzmax);
for j=1:nplots
    axis(hs((jj-1)*nplots+j).h,[-maxXX maxXX 0 maxZZ]);
end;



end;  %figure end

end; %fnmax end
textdat

for jj=2:fnmax
    for j=1:nplots
        (jj-1)*nplots+j
        caxis(hs((jj-1)*nplots+j).h,[0 max(maxval)]);
        colorbar;
        axes('Position',[0 0 1 1],'Visible','off');
        text(.4,.97,textdataDan(1).text);
    end
    bb(jj-1)=getframe(hh(jj));
end


