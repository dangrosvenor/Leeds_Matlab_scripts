%does 2-d plots from multiple directories Q fields 1-9 totalMR=16

clear b xxmax zzmax maxval;
con=input('Enter col no.: ');
cono=con;
nplots=input('Enter the no. of plots:');

scrsz=get(0,'ScreenSize');
posit=[9 10 scrsz(3)/1.01 scrsz(4)/1.1];

owntitle=1;
tit(1).t='Constant Fluxes';
tit(2).t='Sine wave flux variation';
tit(3).t='Double sine wave flux variation';


clear maxims hh;   
jmax=6;

for k=1:fix((nplots-0.1)/6)+1
    hh(k)=figure('name',strcat('2d plots with same scales of column:',int2str(cono)),'Position',posit);
    
    
     if nplots-k*6<0 
         jmax=rem(nplots,6)+(k-1)*6;
     end;
 pp=0;
 
for j=(k-1)*jmax+1:min(jmax,nplots);
    pp=pp+1;
    
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
    
    a=5;
    b=1;
    
    if pp==3
        plotnum=4
    end
    
    hs((k-1)*jmax+j).h=subplot(a,b,plotnum);
    axis equal;
    [ZZ XX]=meshgrid(GridDan(j).Z,GridDan(j).Y1);
    ZZ=ZZ';
    XX=XX';
    xxmax(j)=max(GridDan(j).Y1);
    zzmax(j)=max(GridDan(j).Z);
    
    %for j=1:nplots
    
    
   
    if cono<=9
        pdat=TwoDDan(j).Q(:,:,cono);
        if cono==1;
            titlecb='Vapour MR (kg/kg)';
            pdat(1,:)=0;
        end;
        if cono==2;
            titlecb='Liquid MR (kg/kg)';
            pdat(1,:)=0;
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
        
    end;
    
    if cono==111
        pdat=TwoDDan(j).U;
        titlecb='U (m/s)';
    end;
    if cono==112
        pdat=TwoDDan(j).V;
        titlecb='V (m/s)';
    end;
    if cono==113
        pdat=TwoDDan(j).W;
        titlecb='W (m/s)';
    end;
    if cono==114
        pdat=TwoDDan(j).TH;
        titlecb='Theta (K)';
    end;
    
     if cono==16;
        titlecb='Total LRSGI MR (kg/kg)';
        pdat=TwoDDan(j).Q(:,:,2);
        for kk=3:6;
            pdat=pdat+TwoDDan(j).Q(:,:,kk);
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
        
 
    maxims((k-1)*jmax+j)=max(max(pdat));
        
    
    %pdat(1,1)=maxval;
    pcolor(XX./1000,ZZ./1000,pdat);
    
    VW=TwoDDan(j).V;
    WW=TwoDDan(j).W;
    size(VW);
    %hss=streamslice(XX./1000,ZZ./1000,VW,WW,ans(1)/300);
    %set(hss,'Color','black');
    
    if owntitle==1
        title(tit(j).t);
    else
        title(direcDan(j).dir);
    end
    
    colorbar;
    shading interp;
    
    %set(gca,'XTicklabel','');
    
    
    ylabel( texlabel('Height(km)') );
    
    %if j==jmax
        text(max(xxmax)/1000/1.25,-max(zzmax)/1000/11,titlecb);
        %end
    
    %text(10,-0.5,titlecb,'units','centimeters');
    %end;
 %   b(jj-1)=getframe;
 
end; %subplot end
maxval(k)=max(maxims);
maxXX=max(xxmax)/1000;
maxZZ=max(zzmax)/1000;
for j=1:nplots
    axis(hs((k-1)*jmax+j).h,[-maxXX maxXX 0 14]);
end;



end;  %figure end



for jj=1:k
    for j=1:nplots
        caxis(hs((jj-1)*jmax+j).h,[0 max(maxval)]);
        colorbar('peer',hs((jj-1)+j).h);
        axes('Position',[0 0 1 1],'Visible','off');
        text(.4,.98,textdataDan(1).text);
    end
   % bb(jj-1)=getframe(hh(jj));
end

subplot(5,1,3);
plot(Grid.Y1(2:end-1)./1000,sflux(2).s);
set(gca,'XTicklabel','');
axis tight;
colorbar;
ylabel('SHeat Flux(Wm^-2)');
subplot(5,1,5);
plot(Grid.Y1(2:end-1)./1000,sflux(3).s);
axis tight;
colorbar;

ylabel('SHeat Flux(Wm^-2)');
xlabel(texlabel('Horizontal Displacement(km)'));
%text(max(Grid.Y1)/1000/1.1,315,titlecb);


exdir='field\onesine\nicePICS\';

gcf=hh(1);
set(gcf,'paperpositionmode','auto');
exname=strcat('c:\matlabR12\work\',exdir,'dump25.jpg');	
print(gcf,'-djpeg','-r150',exname);
%close(gcf);

