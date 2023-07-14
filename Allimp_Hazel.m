%fnmin=input('Enter the first diagnostic file number, with directories in direcDan(i).dir: ');
%fnmax=input('Enter the maximum diagnostic file number, with directories in direcDan(i).dir: ');
%files=input('Enter the no. of directories to import %files that want to read in from e.g. enter [1 2 4] for direcDan([1 2 4]).dir:');

fnmin=1; %start file to read in
fnmax=10; %end file to read in
files=[1 2 4]; %files that want to read in from e.g. enter [1 2 4] for direcDan([1 2 4]).dir



for jjfile=fnmin:fnmax  %for1
    
    
    fprintf(1,'jjfile=%d ',jjfile);
    if jjfile<10;
        fn=strcat('00',int2str(jjfile));
    elseif jjfile<100
        fn=strcat('0',int2str(jjfile));
    else
        fn=int2str(jjfile);
    end
    
    fname=['RUN0001.DG0' fn]; %change this to that required by Paul's DIAG program
    
    
    
    for jc=1:length(files)
        idir=files(jc);
        
                
        FileName=[direcDan(idir).dir fname];
        
        clear TwoD ThreeD
        
        
        DIAG_3d_newt;  %%change this to Paul's DIAG program
        
        
        
        if jjfile==fnmin;
            sortheaderDGS; %sorts out the HEADER file to put ALL_Q01 etc. in the right order for using in TimeAv.DAGV
                           %there is prob something in Paul's script that will replace this
            dgstrDan(idir).dg=dgstr;
        end
        
        
        
        
        TwoDDan(idir)=TwoD; 
        
        SerDan(idir).SER=SER;
        
        TimeAvDan(idir).RNDGS=TimeAv.RNDGS;
        TimeAvDan(idir).DGAV=TimeAv.DGAV;
        
        icediags_5thSept_Hazel;
        
        
        
        
    end
    
end
    
    