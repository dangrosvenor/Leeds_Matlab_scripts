function [fields,varargout] = DIAG2_3_3D_LINUX(FileName)
% M -FILE to Import the DIAG file for LES version 2.3
% START import program
fid=fopen(FileName,'rb');jjArgOut = 1;kkArgOut = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just skip the header
%fseek(fid,65.*512.,'bof');
HEADER2='';
while(1)
    HEADER=((fread(fid,512,'char')))';
    HEADER2=cat(2,HEADER2,HEADER);
    if(findstr(HEADER,'END:'))
        break;
    end
end
%[Grid.DG,Grid.SERSTR]=SORTHEADER(HEADER2);
TimeAv.HEADER = HEADER2;
%fseek(fid,-8,'cof');
%HEADER=fread(fid,i,'char=>char');

% APARMS is a float array of length 1024
APARMS=fread(fid,1024,'float=>double');



% Now define the offsets
OFFSET(1)=(46.+APARMS(45));
OFFSET(2)=OFFSET(1) +double(APARMS(17));
OFFSET(3)=OFFSET(1) +double(2.*APARMS(17));
OFFSET(4)=OFFSET(1) +double(3.*APARMS(17));
OFFSET(5)=OFFSET(4) +double(49.);
OFFSET(6)=OFFSET(5) +double(APARMS(17));
IHS=OFFSET(6) +double(15.);
NPARMS=IHS+67;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET A FEW POINTERS
POINT.IPASQP = APARMS(16);
POINT.IRAINP = APARMS(32);
if(POINT.IPASQP == 0 & POINT.IRAINP >= 10)
    POINT.IQV = APARMS(IHS+22);
    POINT.IQL = 2;  % HARD CODED IQL --Think about adding to APARMS ARRAY at some point
    POINT.IQR = APARMS(IHS+23);
    POINT.IQS = APARMS(IHS+24);
    POINT.IQG = APARMS(IHS+25);
    POINT.IQI = APARMS(IHS+26);
    POINT.IQN = APARMS(IHS+27);
    POINT.IQNR = APARMS(IHS+28);
    POINT.IQNS = APARMS(IHS+29);
    POINT.IQNG = APARMS(IHS+30);
    POINT.IQGV = APARMS(IHS+66);
    % FIX FOR THE TRACER POINTERS --NB in future think about adding to aparms array
    if(max([POINT.IQV POINT.IQL POINT.IQR POINT.IQS POINT.IQG ...
                POINT.IQI POINT.IQN POINT.IQNR POINT.IQNS POINT.IQNG POINT.IQGV]) ==  APARMS(17)-1)
        POINT.ITR1 = APARMS(17);
    end
    if(max([POINT.IQV POINT.IQL POINT.IQR POINT.IQS POINT.IQG ...
                POINT.IQI POINT.IQN POINT.IQNR POINT.IQNS POINT.IQNG POINT.IQGV]) ==  APARMS(17)-2)
        POINT.ITR1 = APARMS(17)-1;
        POINT.ITR2 = APARMS(17);
    end
end
% In this version the blocks are blocks of 128
Grid.X1=Read1D(APARMS(1)+2,fid);
Grid.X2=Read1D(APARMS(1)+2,fid);
Grid.Y1=Read1D(APARMS(2)+2,fid);
Grid.Y2=Read1D(APARMS(2)+2,fid);
Grid.Z=Read1D(APARMS(3),fid);
Grid.ZN=Read1D(APARMS(3),fid);
Grid.RHO=Read1D(APARMS(3),fid);
Grid.RHON=Read1D(APARMS(3),fid);

% Choice of IDGx dependent
if(APARMS(OFFSET(5)+1) > 0)
    Grid.UBAR=Read1D(APARMS(3),fid);
end
if(APARMS(OFFSET(5)+2) > 0)
    Grid.VBAR=Read1D(APARMS(3),fid);
end
if(APARMS(OFFSET(5)+5) > 0)
    Grid.OLTHBAR=Read1D(APARMS(3),fid);
    Grid.THINIT=Read1D(APARMS(3),fid);
    Grid.THREF=Read1D(APARMS(3),fid);
end
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) > 0)
        Grid.OLQBAR(:,IQ)=Read1D(APARMS(3),fid);
    end
end
if(APARMS(OFFSET(5)+6) > 0)
    Grid.PREFN=Read1D(APARMS(3),fid);
end

% Now XXX
if(APARMS(35) >= 1 & APARMS(36) ~= 0)
    % XXX_A
    for IL=1:APARMS(35)
        Grid.XXX_A(:,IL)=Read1D(APARMS(3),fid);
    end
    
    % XXX_W
    for IL=1:APARMS(35)
        Grid.XXX_W(:,IL)=Read1D(APARMS(3),fid);
    end
        
    % XXX_TH
    for IL=1:APARMS(35)
        Grid.XXX_TH(:,IL)=Read1D(APARMS(3),fid);
    end
        
    % XXX_QYY
    for IQ=1:APARMS(17)
        for IL=1:APARMS(35)
            Grid.XXX_QYY(:,IQ,IL)=Read1D(APARMS(3),fid);
        end
    end
end

% Import the time-averaged diagnostics
if(APARMS(OFFSET(6)+8) ==1)
    TimeAv.RNDGS=Read1D(APARMS(OFFSET(6)+9),fid);
    TimeAv.DGAV=Read1D2(APARMS(3),APARMS(OFFSET(6)+9),fid);
    %TimeAv.DGAV=Read1D(APARMS(OFFSET(6)+9).*APARMS(3),fid);
    varargout{jjArgOut} = TimeAv;jjArgOut = jjArgOut + 1;
    fields{kkArgOut} = 'TimeAv';kkArgOut = kkArgOut + 1;
end

if(APARMS(OFFSET(6)+14) ==1 & APARMS(OFFSET(6)+13) >=1)
    %Grid.DGSPAV=Read1D2(APARMS(OFFSET(6)+13),floor(JJTOTP./2+1),fid);
    Grid.DGSPAV=Read1D(APARMS(OFFSET(6)+13).*(APARMS(44)./2+1),fid);
    
end
% Time Series
if(APARMS(28) == 1)
    for NSER=1:APARMS(34)
        SER(:,NSER)=Read1D(APARMS(OFFSET(6)+10),fid);
    end
    varargout{jjArgOut} = SER;jjArgOut = jjArgOut + 1;
    fields{kkArgOut} = 'SER';kkArgOut = kkArgOut + 1;
end
if(exist('SER','var') & exist('SER','var'))
    Grid.TIME(:,1)=squeeze(SER(end,1));
end
if(exist('Grid','var'))
    varargout{jjArgOut} = Grid;jjArgOut = jjArgOut + 1;
    fields{kkArgOut} = 'Grid';kkArgOut = kkArgOut + 1;
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
    TwoD.U=Read2D(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+2) == 2)
    TwoD.V=Read2D(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+3) == 2)
    TwoD.W=Read2D(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+5) == 2)
    TwoD.TH1=Read2D(APARMS(2)+2,APARMS(3),fid);
end
% 2-D Q variables
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) == 2)
        TwoD.Q(:,:,IQ)=Read2D(APARMS(2)+2,APARMS(3),fid);
    end
end
% RADIATION HERE ***********************
if(APARMS(37)>=1)
    TwoD.sthrad=Read2D(APARMS(2)+2,APARMS(3),fid);
    if(APARMS(38)==1)
        TwoD.twdrad=Read2D(APARMS(2)+2,APARMS(3),fid);
    end
    if(APARMS(39)==1)
        TwoD.swdrad=Read2D(APARMS(2)+2,APARMS(3),fid);
    end
end
% more code here
    if(APARMS(37)>=2)
        if(APARMS(38)==1)
            TwoD.olr=Read2D(APARMS(1)+2,APARMS(2)+2,fid);
            TwoD.dslw=Read2D(APARMS(1)+2,APARMS(2)+2,fid);
            TwoD.uslw=Read2D(APARMS(1)+2,APARMS(2)+2,fid);
        end
        if(APARMS(39)==1)
            TwoD.dtsw=Read2D(APARMS(1)+2,APARMS(2)+2,fid);
            TwoD.utsw=Read2D(APARMS(1)+2,APARMS(2)+2,fid);
            TwoD.dssw=Read2D(APARMS(1)+2,APARMS(2)+2,fid);
            TwoD.ussw=Read2D(APARMS(1)+2,APARMS(2)+2,fid);
        end
        if(APARMS(37)>=3)
            if(APARMS(38)==1)
                TwoD.swfxup=Read2D(APARMS(2)+2,APARMS(3),fid);
                TwoD.swfxdt=Read2D(APARMS(2)+2,APARMS(3),fid);
                TwoD.swfxnt=Read2D(APARMS(2)+2,APARMS(3),fid);
            end
            
            if(APARMS(39)==1)
                TwoD.twfxup=Read2D(APARMS(2)+2,APARMS(3),fid);
                TwoD.twfxnt=Read2D(APARMS(2)+2,APARMS(3),fid);
            end
        end
    end
            
if(APARMS(OFFSET(5)+6) == 2)
    TwoD.P=Read2D(APARMS(2)+2,APARMS(3),fid);
    RHOref=repmat(Grid.RHON,[1 length(Grid.Y1)]);
    Pref=repmat(Grid.PREFN,[1 length(Grid.Y1)]);
    TwoD.PP=TwoD.P.*RHOref + Pref;
end
if(APARMS(OFFSET(5)+7) == 2)
    TwoD.PV=Read2D(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(40) == 2)
    TwoD.TH2=Read2D(APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(41) == 2)
    TwoD.THV=Read2D(APARMS(2)+2,APARMS(3),fid);
end
if(exist('TwoD','var'))
    varargout{jjArgOut} = TwoD;jjArgOut = jjArgOut + 1;
    fields{kkArgOut} = 'TwoD';kkArgOut = kkArgOut + 1;
end

% 3-d
if(APARMS(OFFSET(5)+1)==3)
    ThreeD.U=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+2)==3)
    ThreeD.V=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+3)==3)
    ThreeD.W=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+5)==3)
    ThreeD.TH1=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
% Q
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) == 3)
        ThreeD.Q(:,:,:,IQ)=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
    end
end

if(APARMS(OFFSET(5)+6) == 3)%GP
    ThreeD.P=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(OFFSET(5)+7) == 3)%PV
    ThreeD.PV=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(40)==3)%theta
    ThreeD.TH2=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
if(APARMS(41)==3)%THV
    ThreeD.THV=Read3D(APARMS(1)+2,APARMS(2)+2,APARMS(3),fid);
end
if(exist('ThreeD','var'))
    varargout{jjArgOut} = ThreeD;jjArgOut = jjArgOut + 1;
    fields{kkArgOut} = 'ThreeD';kkArgOut = kkArgOut + 1;
end
%NMET BIT HERE
%fseek(fid,-4.*128,'eof');
    TEMP=Read1D(APARMS(3),fid);

c=ftell(fid);
fseek(fid,0,'eof');
d=ftell(fid);
(d-c)./4;
%clear IHS;
%clear IQ;
%clear IL;
%clear NSER;
%clear NPARMS;

fclose(fid);
%clear fid;
varargout{jjArgOut} = POINT;jjArgOut = jjArgOut + 1;
fields{kkArgOut} = 'POINT';;kkArgOut = kkArgOut + 1;

if(nargout ~= jjArgOut)
    error(sprintf('Must have %d outputs',jjArgOut));
end
