%Run Diag_3d_newt with idontclose=1 to stop the import just before the Q-field1 are due to be read in

clear TwoD
    
fnall=input('Enter the Diagnostic file number, with directories in direcDan(i).dir: ');
ndir=input('Enter the no. of directories to import :');

dirs=ndir(2:end);
ndir=ndir(1);

sorthead=1;
    

size(fnall);
fnsiz=ans(2);

for jj=1:ndir;
    
    if length(dirs)>=1 %if want to import only certain files then enter no. files and then file numbers to be imported
        idir=dirs(jj); %e.g. enter [3 5 6 7] to enter files 5,6&7
    else
        idir=jj;
    end
    
    if fnsiz>1
        fn=fnall(idir);
    else
        fn=fnall;
    end
    
    if fn<10;
        fn=strcat('00',int2str(fn));
    elseif fn<100
        fn=strcat('0',int2str(fn));
    else
        fn=int2str(fn);
    end
    
	fname=['RUN0001.DG0' fn];    
    FileName=[direcDan(idir).dir fname];   
    
    fprintf(1,'\n Loading file...');
    Diag_3d_newt;
    try
    GridDan(jj)=Grid;
	catch
        fprintf(1,'\n**** Error with "GridDan(jj)=Grid;"****\n');
    end
    fprintf(1,'\n3d slice getter running.....');
    f=1e6*28.97/18;
    imax=APARMS(1)+2;
    jmax=NJ+2;
    kmax=APARMS(3);
    
    
    
    ixpos=152;
    ixpos=2;
    
    %ixpos=2;
    %iypos=2;
    
    
    
    
    for iq=1:9   
        
        for k=1:kmax             
            
            %read in new value for iq
            X=fread(fid,[imax.*jmax],'float=>double');           
            X=reshape(X,[imax jmax]);
           % X=X(2:end-1,:);
            
            TwoDDan(idir).Q(k,1:size(X,1),iq)=X(:,ixpos);
          %  	TwoD.Q(k,:,iq)=X(:,ixpos)';                                  
        end
        
    end
    
end
  
%TwoDDan(idir)=TwoD;

%GridDan(idir)=Grid;
                    
		
            
    fprintf(1,'\ndone');
%    break
    
