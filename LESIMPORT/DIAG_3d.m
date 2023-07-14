% M -FILE to Import the DIAG file for LES version 2.3
% START import program
clear SER TwoD ThreeD oneD TimeAv surfdg HEADER2;

bezier=0;
newton=0;
if bezier==1
    fid=fopen(FileName,'rb','s');
    linux=1;
elseif newton==1
    fid=fopen(FileName,'rb','a');
    linux=2;
else
    fid=fopen(FileName,'rb');
    linux=0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just skip the header
%fseek(fid,18.*512.*4,'bof');

if linux==0

hsize=512;
ind=0;
HEADER=((fread(fid,hsize,'2048*char')))';

if linux==1
    ha=1024;
else
    ha=2048;
end


fseek(fid,ha-hsize,'cof');
   
%char(HEADER)
textdat=char(HEADER(35:48));
ihead=1;
HEADER2((ihead-1)*hsize+1:ihead*hsize)=HEADER;
ihead=2;
while(1) & ihead<500
    HEADER=((fread(fid,hsize,'2048*char')))';
    fseek(fid,ha-hsize,'cof');  
    icomma=1;
    b=findstr('I3&3',HEADER);
    if(b)
        'FOUND AN I3&3 - COULD BE ERRORS IN SORTING OF HEADER'
        %if icomma==1; HEADER(b(1))=','; end%add comma - should only do this in certain cases where I3&3 has written over the original comma
        aa=findstr(',',HEADER);
        ad=aa(2:end)-aa(1:end-1);
        auni=unique(ad);
        if length(auni)>1
            HEADER(b(1))=','; %add comma if I3&3 has written over the previous comma and caused a diff in length between commas
        end                     %**** CHECK to see whether has written over name of anything useful *****
    end
    
    
    %char(HEADER)
    if(findstr(HEADER,'DGEND')) %will have to change depending on location of split
            HEADER2((ihead-1)*hsize+1:ihead*hsize)=HEADER;
        break;
        
    end
    
    
    
    HEADER2((ihead-1)*hsize+1:ihead*hsize)=HEADER;
    ihead=ihead+1;
end

dg=findstr(HEADER,'DGEND:');

if findstr(HEADER,'DGEND') & findstr(HEADER,'DGEND:') 
else
    fseek(fid,2048,'cof'); %if DGEND: is split over successive records will have to change depending on location of split
end


else

ihead=0;
HEADER2='';
while(1) & ihead<2000
    ihead=ihead+1;
    HEADER=((fread(fid,512,'char')))';
    HEADER2=cat(2,HEADER2,HEADER);
    if(findstr(HEADER,'END:'))
        break;
    end
end

end

textdat=char(HEADER2(35:48));

if ihead==2000
    'SOMETHING WRONG WITH FINDING "END"'
    break
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

%test=fread(fid,512,'2048*char');

if newton==1
fseek(fid,1024+128*4,'cof');
end

%APARMS=fread(fid,1024,'1024*float=>double');
APARMS=fread(fid,1024,'float=>double');


% Now define the offsets
OFFSET(1)=(45.+APARMS(45)); %=ISTART1
OFFSET(2)=OFFSET(1) +double(APARMS(17));
OFFSET(3)=OFFSET(1) +double(2.*APARMS(17));
OFFSET(4)=OFFSET(1) +double(3.*APARMS(17));
OFFSET(5)=OFFSET(4) +double(49.);   %2nd ISTART2=ISTART+49+3*NQP
OFFSET(6)=OFFSET(5) +double(APARMS(17));
IHS=OFFSET(6) +double(15.);
NPARMS=IHS+67;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if newton==1
fseek(fid,1024*3*4,'cof');
end


% if linux==1
% 	fseek(fid,4*1024,'cof');
% end

part=128;
% In this version the blocks are blocks of 128
[Grid.X1 part]=Read1D_vrad(APARMS(1)+2,fid,linux,part);
[Grid.X2 part]=Read1D_vrad(APARMS(1)+2,fid,linux,part);


if APARMS(2)>1500  %needed to add for JJP=1536 run. not sure where threshold lies
%fseek(fid,1024*4,'cof');
end

Grid.Y1=Read1D_vrad(APARMS(2)+2,fid,linux);
Grid.Y2=Read1D_vrad(APARMS(2)+2,fid,linux);
Grid.Z=Read1D_vrad(APARMS(3),fid,linux);
Grid.ZN=Read1D_vrad(APARMS(3),fid,linux);



Grid.RHO=Read1D_vrad(APARMS(3),fid,linux);


Grid.RHON=Read1D_vrad(APARMS(3),fid,linux);


% Choice of IDGx dependent
if(APARMS(OFFSET(5)+1) > 0)
    Grid.UBAR=Read1D_vrad(APARMS(3),fid,linux);
end
if(APARMS(OFFSET(5)+2) > 0)
    Grid.VBAR=Read1D_vrad(APARMS(3),fid,linux);
end
if(APARMS(OFFSET(5)+5) > 0)
    Grid.OLTHBAR=Read1D_vrad(APARMS(3),fid,linux);
    Grid.THINIT=Read1D_vrad(APARMS(3),fid,linux);
    Grid.THREF=Read1D_vrad(APARMS(3),fid,linux);
end
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) > 0)
        Grid.OLQBAR(1:APARMS(3),IQ)=Read1D_vrad(APARMS(3),fid,linux);
    end
end
if(APARMS(OFFSET(5)+6) > 0)
    Grid.PREFN=Read1D_vrad(APARMS(3),fid,linux);
end

IPRAB=NPARMS;
NITEST=APARMS(IPRAB+1);
NJTEST=APARMS(IPRAB+2);
NKTEST=APARMS(IPRAB+3);

NPRAB=IPRAB+3+NITEST+NJTEST+NKTEST+1+2+5; %postition in APARMS of end of radiation flags and start of flux flags

imaxrad=APARMS(NPRAB+8);
itrtop=APARMS(NPRAB+9);
imaxpr=APARMS(NPRAB+10);
icape=APARMS(NPRAB+12);
iaero=APARMS(NPRAB+20);

itrtop=1;

 if APARMS(NPRAB+1)==1 % read in surface varying fluxes

     SENflux=Read1D_vrad(APARMS(2),fid,linux).*Grid.RHO(1).*1005; %convert to W/m2
     LATflux=Read1D_vrad(APARMS(2),fid,linux).*Grid.RHO(1).*2.501e6;
 end
 
 if APARMS(NPRAB+2)==1 % read in surface varying roughness lengths

     Z0var=Read1D_vrad(APARMS(2),fid,linux);
     Z0THvar=Read1D_vrad(APARMS(2),fid,linux);
 end



% Now XXX
if(APARMS(35) >= 1 & APARMS(36) ~= 0)
    % XXX_A
    for IL=1:APARMS(35)
        Grid.XXX_A(1:APARMS(3),IL)=Read1D_vrad(APARMS(3),fid,linux);
    end
    
    % XXX_W
    for IL=1:APARMS(35)
        Grid.XXX_W(1:APARMS(3),IL)=Read1D_vrad(APARMS(3),fid,linux);
    end
        
    % XXX_TH
    for IL=1:APARMS(35)
        Grid.XXX_TH(1:APARMS(3),IL)=Read1D_vrad(APARMS(3),fid,linux);
    end
        
    % XXX_QYY
    for IQ=1:APARMS(17)
        for IL=1:APARMS(35)
            Grid.XXX_QYY(1:APARMS(3),IQ,IL)=Read1D_vrad(APARMS(3),fid,linux);
        end
    end
end
% Import the time-averaged diagnostics
if(APARMS(OFFSET(6)+8) ==1)
    TimeAv.RNDGS=Read1D_vrad(APARMS(OFFSET(6)+9),fid,linux);
    TimeAv.DGAV=Read1D2_vrad(APARMS(3),APARMS(OFFSET(6)+9),fid,linux);
    %TimeAv.DGAV=Read1D(APARMS(OFFSET(6)+9).*APARMS(3),fid);
end


if(APARMS(NPRAB+11) >0)
    for k=1:APARMS((NPRAB+11))
        TimeAv.surav(:,:,k)=Read2D_vrad(APARMS(2),APARMS(1),fid,linux);
    end
end



if(APARMS(OFFSET(6)+14) ==1 & APARMS(OFFSET(6)+13) >=1)
    %Grid.DGSPAV=Read1D2(APARMS(OFFSET(6)+13),floor(JJTOTP./2+1),fid,linux);
    Grid.DGSPAV=Read1D_vrad(APARMS(OFFSET(6)+13).*(APARMS(44)./2+1),fid,linux);
end
% Time Series
if(APARMS(28) == 1)
   for NSER=1:APARMS(34)
        SER(:,NSER)=Read1D_vrad(APARMS(OFFSET(6)+10),fid,linux);
end
end

if (APARMS(37)==0) %if no radiation
    SER(:,38:APARMS(34)+7)=SER(:,31:APARMS(34));  %move timeseries down so that process rates are in the places listed on sheet
    SER(:,31:37)=0;
end

nqp=APARMS(17);
diff=nqp-9;

endno=2*imaxrad+itrtop+2*imaxpr+2*icape+2*iaero;
%col 108=maxrad 109=maxrad height 110=max tracer height
%col 111=max precipitation rate 112=J location of max pp rate
%col 113=average cape, 114=av CIN
%115=max supersat (ev/es) (Nenes scheme), 116=height of max

nclear=500;
if (nqp>9)
    SER(:,nclear:nclear-1+diff)=SER(:,38+nqp-diff:38+nqp-1); %put extra QxxMAX's for Q's above 9 into nclear+ columns
    SER(:,nclear+diff:nclear-1+(2*diff))=SER(:,38+(2*nqp)-diff:38+(2*nqp)-1); %extra HQxxMAX's
    
               %SER(:,nclear+(2*diff)):nclear-1+(3*diff)))=SER(:,38+(2*nqp)+3-diff+1:38+(2*nqp)+3);
    
    SER(:,47:55)=SER(:,38+nqp:38+nqp+8);
    SER(:,56:107+endno)=SER(:,38+(2*nqp):38+(2*nqp)+51+endno);
    
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
    ThreeD=0;
    TwoD.U=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
end
if(APARMS(OFFSET(5)+2) == 2)
    TwoD.V=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
end
if(APARMS(OFFSET(5)+3) == 2)
    TwoD.W=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
end
if(APARMS(OFFSET(5)+5) == 2)
    TwoD.TH1=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
end
% 2-D Q variables
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) == 2)
        TwoD.Q(:,:,IQ)=Read2D_vradNOSKIP(APARMS(2)+2,APARMS(3),fid,linux);
    end
end









% RADIATION HERE ***********************
if(APARMS(37)>=1)
    TwoD.sthrad=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
    if(APARMS(38)==1)
        TwoD.twdrad=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
    end
    if(APARMS(39)==1)
        TwoD.swdrad=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
    end
end
% more code here
   % if(APARMS(37)>=2)
   %     if(APARMS(38)==1)

   
if(APARMS(OFFSET(5)+6) == 2)
    TwoD.P=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
    RHOref=repmat(Grid.RHON,[1 length(Grid.Y1)]);
    Pref=repmat(Grid.PREFN,[1 length(Grid.Y1)]);
    TwoD.PP=TwoD.P.*RHOref + Pref;
else
    TwoD.P=0;
    TwoD.PP=0;
end
if(APARMS(OFFSET(5)+7) == 2)
    TwoD.PV=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
else
    TwoD.PV=0;
end
if(APARMS(40) == 2)
    TwoD.TH2=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
else
    TwoD.TH2=0;
end
if(APARMS(41) == 2)
    TwoD.THV=Read2D_vrad(APARMS(2)+2,APARMS(3),fid,linux);
else
    TwoD.THV=0;
end


if(APARMS(OFFSET(5)+6) == 2)

%     Pref=repmat(Grid.PREFN./Grid.RHON,[1 length(Grid.Y1)]);
%     Pdash=Pref+TwoD.P; %pressure perturbation added to reference P/rho
%     TwoD.PP=(Pdash.*28.97e-3/8.31./TwoD.TH2).^(1/0.286)*1e5; %actual pressure

%P field is P'/rhoREF!!!!and NOT (p/rho)'

end



%NMET BIT HERE
%fseek(fid,-4.*128,'eof');
    %TEMP=Read1D_vrad(APARMS(3),fid);
    
if(APARMS(NPRAB+13)==0) %igaloffp=0 - Galilean transform
    surfdiag.precgal=Read2D_vrad(APARMS(NPRAB+16),APARMS(NPRAB+15),fid,linux); 
    surfdiag.pmaxgal=Read2D_vrad(APARMS(NPRAB+16),APARMS(NPRAB+15),fid,linux); 
    surfdiag.instantgal=Read2D_vrad(APARMS(NPRAB+16),APARMS(NPRAB+15),fid,linux); 
end
    

if(APARMS(NPRAB+14)==1)
    for imetp=1:APARMS(45)
	    if(APARMS(45+imetp)==2)
            surfdiag.prec(1:APARMS(1),1:APARMS(2),imetp)=Read2D_vrad(APARMS(2),APARMS(1),fid,linux);    %total acc precip over ALL of runtime
	    end
    end
    surfdiag.pmax=Read2D_vrad(APARMS(2),APARMS(1),fid,linux); %max prate over diag time
    surfdiag.class=Read2D_vrad(APARMS(2),APARMS(1),fid,linux); 
    surfdiag.citop=Read2D_vrad(APARMS(2),APARMS(1),fid,linux); 
    surfdiag.lwp=Read2D_vrad(APARMS(2),APARMS(1),fid,linux);
    surfdiag.iwp=Read2D_vrad(APARMS(2),APARMS(1),fid,linux);
    surfdiag.instant=Read2D_vrad(APARMS(2),APARMS(1),fid,linux);
end


xlen=APARMS(2)+2;


if(APARMS(OFFSET(5)+1) == 3)
     for i=1:xlen^2
         fseek(fid,8,'cof');
         fseek(fid,4,'cof');
    end
    
    for j=1:xlen
            fseek(fid,8,'cof');
            fseek(fid,8,'cof');
            fseek(fid,8,'cof');
            fseek(fid,8,'cof');
            fseek(fid,8,'cof');
            fseek(fid,8,'cof');
    end
    
    ThreeD.U=Read3D2(APARMS(2)+2,APARMS(2)+2,APARMS(3),fid,linux);
end



if(APARMS(OFFSET(5)+2) == 3)
    ThreeD.V=Read3D2(APARMS(2)+2,APARMS(2)+2,APARMS(3),fid,linux);
end

if(APARMS(OFFSET(5)+3) == 3)
    ThreeD.W=Read3D2(APARMS(2)+2,APARMS(2)+2,APARMS(3),fid,linux);
end
if(APARMS(OFFSET(5)+5) == 3)
    ThreeD.TH1=Read3D2(APARMS(2)+2,APARMS(2)+2,APARMS(3),fid,linux);
end
% 2-D Q variables
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) == 3)
        ThreeD.Q(:,:,:,IQ)=Read3D2(APARMS(2)+2,APARMS(2)+2,APARMS(3),fid,linux);
    end
end

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

