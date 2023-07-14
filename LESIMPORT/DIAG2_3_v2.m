% M -FILE to Import the DIAG file for LES version 2.3
% START import program
clear SER TwoD;
fid=fopen(FileName,'rb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just skip the header
%fseek(fid,18.*512.*4,'bof');

hsize=512;
ind=0;
HEADER=((fread(fid,512,'2048*char')))';
fseek(fid,2048-512,'cof');
%char(HEADER)
textdat=char(HEADER(35:48));
ihead=1;
HEADER2((ihead-1)*hsize+1:ihead*hsize)=HEADER;
ihead=2;
while(1)
    HEADER=((fread(fid,512,'2048*char')))';
    fseek(fid,2048-512,'cof');
    %char(HEADER)
    if(findstr(HEADER,'END:'))
            HEADER2((ihead-1)*hsize+1:ihead*hsize)=HEADER;
        break;
        
    end
    HEADER2((ihead-1)*hsize+1:ihead*hsize)=HEADER;
    ihead=ihead+1;
end
%i=1;charc=1;
%while(charc~=256)
%   charc=fread(fid,1,'float=>double');
%    i=i+1;
%end
%i;
%fseek(fid,-8,'cof');
%HEADER=fread(fid,i,'char=>char');

% APARMS is a float array of length 1024
APARMS=fread(fid,1024,'1024*float=>double');

% Now define the offsets
OFFSET(1)=(45.+APARMS(45));
OFFSET(2)=OFFSET(1) +double(APARMS(17));
OFFSET(3)=OFFSET(1) +double(2.*APARMS(17));
OFFSET(4)=OFFSET(1) +double(3.*APARMS(17));
OFFSET(5)=OFFSET(4) +double(49.);
OFFSET(6)=OFFSET(5) +double(APARMS(17));
IHS=OFFSET(6) +double(15.);
NPARMS=IHS+67;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% In this version the blocks are blocks of 128
Grid.X1=Read1D_vrad(APARMS(1)+2,fid);
Grid.X2=Read1D_vrad(APARMS(1)+2,fid);


if APARMS(2)>1500  %needed to add for JJP=1536 run. not sure where threshold lies
fseek(fid,1024*4,'cof');
end

Grid.Y1=Read1D_vrad(APARMS(2)+2,fid);
Grid.Y2=Read1D_vrad(APARMS(2)+2,fid);
Grid.Z=Read1D_vrad(APARMS(3),fid);
Grid.ZN=Read1D_vrad(APARMS(3),fid);



Grid.RHO=Read1D_vrad(APARMS(3),fid);


Grid.RHON=Read1D_vrad(APARMS(3),fid);


% Choice of IDGx dependent
if(APARMS(OFFSET(5)+1) > 0)
    Grid.UBAR=Read1D_vrad(APARMS(3),fid);
end
if(APARMS(OFFSET(5)+2) > 0)
    Grid.VBAR=Read1D_vrad(APARMS(3),fid);
end
if(APARMS(OFFSET(5)+5) > 0)
    Grid.OLTHBAR=Read1D_vrad(APARMS(3),fid);
    Grid.THINIT=Read1D_vrad(APARMS(3),fid);
    Grid.THREF=Read1D_vrad(APARMS(3),fid);
end
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) > 0)
        Grid.OLQBAR(:,IQ)=Read1D_vrad(APARMS(3),fid);
    end
end
if(APARMS(OFFSET(5)+6) > 0)
    Grid.PREFN=Read1D_vrad(APARMS(3),fid);
end

% Now XXX
if(APARMS(35) >= 1 & APARMS(36) ~= 0)
    % XXX_A
    for IL=1:APARMS(35)
        Grid.XXX_A(:,IL)=Read1D_vrad(APARMS(3),fid);
    end
    
    % XXX_W
    for IL=1:APARMS(35)
        Grid.XXX_W(:,IL)=Read1D_vrad(APARMS(3),fid);
    end
        
    % XXX_TH
    for IL=1:APARMS(35)
        Grid.XXX_TH(:,IL)=Read1D_vrad(APARMS(3),fid);
    end
        
    % XXX_QYY
    for IQ=1:APARMS(17)
        for IL=1:APARMS(35)
            Grid.XXX_QYY(:,IQ,IL)=Read1D_vrad(APARMS(3),fid);
        end
    end
end
% Import the time-averaged diagnostics
if(APARMS(OFFSET(6)+8) ==1)
    TimeAv.RNDGS=Read1D_vrad(APARMS(OFFSET(6)+9),fid);
    TimeAv.DGAV=Read1D2_vrad(APARMS(3),APARMS(OFFSET(6)+9),fid);
    %TimeAv.DGAV=Read1D(APARMS(OFFSET(6)+9).*APARMS(3),fid);
end

if(APARMS(OFFSET(6)+14) ==1 & APARMS(OFFSET(6)+13) >=1)
    %Grid.DGSPAV=Read1D2(APARMS(OFFSET(6)+13),floor(JJTOTP./2+1),fid);
    Grid.DGSPAV=Read1D_vrad(APARMS(OFFSET(6)+13).*(APARMS(44)./2+1),fid);
end
% Time Series
if(APARMS(28) == 1)
   for NSER=1:APARMS(34)
        SER(:,NSER)=Read1D_vrad(APARMS(OFFSET(6)+10),fid);
end
end

if (APARMS(37)==0) %if no radiation
    SER(:,38:APARMS(34)+7)=SER(:,31:APARMS(34));  %move timeseries down so that process rates are in the places listed on sheet
    SER(:,31:37)=0;
end

nqp=APARMS(17);
diff=nqp-9;

if (nqp>9)
    SER(:,200:200-1+diff)=SER(:,38+nqp-diff:38+nqp-1); %put extra QxxMAX's for Q's above 9 into 200+ columns
    SER(:,200+diff:200-1+(2*diff))=SER(:,38+(2*nqp)-diff:38+(2*nqp)-1); %extra HQxxMAX's
    %SER(:,200+(2*diff)):200-1+(3*diff)))=SER(:,38+(2*nqp)+3-diff+1:38+(2*nqp)+3);
    
    SER(:,47:55)=SER(:,38+nqp:38+nqp+8);
    SER(:,56:107)=SER(:,38+(2*nqp):38+(2*nqp)+51);
    
end

if (nqp<9)
    SER(:,56:107)=SER(:,38+(2*nqp):38+(2*nqp)+51); %move timeseries so process rate columns are correct
    SER(:,47:47+nqp-1)=SER(:,38+nqp:38+(2*nqp)-1); %shift hqxxMAX's up so go from 47
end
    
    
% Time Series
%if(APARMS(28) == 1)
%   for NSER=1:APARMS(OFFSET(6)+10)
%        SER(:,NSER)=Read1D(APARMS(27),fid);
%    end
%end
%fseek(fid,-4.*(258.*2+122+258),'cof');
% 2-D arrays only valid for 1 slice ???
if(APARMS(OFFSET(5)+1) == 2)
    TwoD.U=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+2) == 2)
    TwoD.V=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+3) == 2)
    TwoD.W=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+5) == 2)
    TwoD.TH1=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
end
% 2-D Q variables
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) == 2)
        TwoD.Q(:,:,IQ)=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
    end
end
% RADIATION HERE ***********************
if(APARMS(37)>=1)
    TwoD.sthrad=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
    if(APARMS(38)==1)
        TwoD.twdrad=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
    end
    if(APARMS(39)==1)
        TwoD.swdrad=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
    end
end
% more code here
   % if(APARMS(37)>=2)
   %     if(APARMS(38)==1)
            
    

if(APARMS(OFFSET(5)+7) == 2)
    TwoD.PV=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(40) == 2)
    TwoD.TH2=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(41) == 2)
    TwoD.THV=Read2D_vrad(APARMS(2)+2,APARMS(3),fid);
end

%NMET BIT HERE
%fseek(fid,-4.*128,'eof');
    TEMP=Read1D_vrad(APARMS(3),fid);

%c=ftell(fid)
fseek(fid,0,'eof');
%d=ftell(fid)
%(d-c)./4
%IHS;
clear IQ;
clear IL;
clear NSER;
clear NPARMS;

fclose(fid);
clear fid;

