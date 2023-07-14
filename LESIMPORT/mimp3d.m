%imports n files directly from .DG00xx files from n directories in direcDan(i).dir
%put the directory paths to RUN* files in direcDan(i).dir where i represents different runs
%enter number of directories to import from if are contained in direcDan(1:n).dir or enter e.g. [3 1 3 5] to load in 3 files (first number)
%for i=1,3&5 in direcDan(i).dir (next n numbers).
%need to set up an array called "machine" to say what type of machine the runs were done on - see DIAG_newt_3d.m for details
%loads in 2-d fields & time averages
%puts 2-d fields into TwoDDan(i) structure - contains U,V & W wind fields, TH1 - potential temperature perturbation (above THREF),
%Q - the scalar fields of vapourMR,liquidMR,rainMR,snowMR,graupelMR,iceMR,iceNC,graupelNC and snowNC MR in kg/kg and NC in #/kg,
%sthrad,twdrad&swdrad - radiative heating
%P=pressure perturbation = (p/rho)' PP=actual pressure (Pa)
%PV= potential vorticity
%also loads time averaged profiles into TimeAvDan(i)
%structure contains DGAV, which contains all the profiles in nzpoints x nprofiles array
%can get column numbers for specific profiles using e.g. col=getDGcol('ALL_Q01',dgstrDan(i).dg); which just looks for string in sorted header
%so then profile of time averaged vapour in this case will be in TimeAvDan(i).DGAV(:,col) - see LEM code for diag descriptions - SUBROUTINE RESDGS

clear SER TwoD SENflux LATflux Z0var Z0THvar dir
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
        fn=fnall(jj);
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

if loadselect==26
	fname=['RUN0003.DG0' fn];
else
	fname=['RUN0001.DG0' fn];
end

    if strcmp(direcDan(idir).dir(2),':')==1
        FileName=[direcDan(idir).dir fname];
    else
        FileName=['c:\cygwin\home\user\' direcDan(idir).dir fname];
    end
    
    if loadselect==26
		 [fields,TimeAv,SER,Grid,TwoD,pointers]=DIAG2_3_3D_LINUX(FileName);
	else
		 DIAG_3d_newt;
	end
    
    pnames; %gives a structure detailing the column numbers of the microphysical process codes for timeseries 
   
    
    if ~exist('sorthead');sorthead=1;end
    if sorthead==1;
        sortheaderDGS;
        dgstrDan(idir).dg=dgstr;
	end
           
    try 
        size(TwoD.Q);
        TwoDDan(idir)=TwoD;
        TwoDDan(idir).Q(:,1,1)=0;
    catch
    end


    SerDan(idir).SER=SER;
    textdataDan(idir).text=textdat;
    
    try 
        TimeAvDan(idir).surav=TimeAv.surav;
    catch
    end   
    TimeAvDan(idir).RNDGS=TimeAv.RNDGS;
    TimeAvDan(idir).DGAV=TimeAv.DGAV;
    
    if length(Grid.X1)>3 & exist('ThreeD')
        ThreeDDan(idir)=ThreeD; %avoiding this to save memory
    end
    
    if ( exist('SENflux') & exist('LATflux') ) 
        sflux(idir).s=SENflux;
        lflux(idir).l=LATflux;
    end
    if ( exist('Z0var') & exist('Z0THvar') ) 
        Z0(idir).z=Z0var;
        Z0TH(idir).z=Z0THvar;
    end
    if exist('surfdiag')
        surfdg(idir)=surfdiag;
    end
    
    try 
        Grid.t; 
    catch
        Grid.t=19.75+SER(end,1)/3600; %set up a dummy Grid.t of the time of the dump if not present already
    end
    
    GridDan(idir)=Grid;
    
end;

sorthead=0;
%[DG,SERSTR]=SORTHEADER(HEADER2);