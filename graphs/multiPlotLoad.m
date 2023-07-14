dirs=1;
nfiles=1;
direcDan(1).dir='g:runs/dmi1715_5ppmv';
exdir='g:runs/dmi1715_5ppmv/results/icesatmr+icenc/';

fnmin=1;
fnmax=88;

for jj=fnmin:fnmax
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
        
        z=Grid.Z;
       
       pressure(i).p(:,:,jj)=TwoDDan(i).PP;
       icenc(i).i(:,:,jj)=TwoDDan(i).Q(:,:,7);
       icemr(i).i(:,:,jj)=TwoDDan(i).Q(:,:,6);
       
       tt=jj;
       plotTimeHeightVap2;
   
       exname=[exdir num2str(jj) '-IceSatMR+iceNC_TIME=' num2str(SER(end,1)/3600) '.emf'];
       set(gcf,'paperpositionmode','auto');
       print(gcf,'-dmeta',exname);
       close(gcf);
       
       
        
end;
end
  