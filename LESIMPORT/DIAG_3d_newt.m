% M -FILE to Import the DIAG file for LES version 2.3
% START import program
clear SER TwoD ThreeD oneD TimeAv surfdg HEADER2;

idontclose=0;   %flag to stop reading in just before Q fields for seperate retrieval of Q
iovride=0;      %flag to set own nreads
%for vertical slice retrieval set both of these to one

if ~exist('ijustheader'); ijustheader=0; end

if iovride==0
	in=1;
	nreads(in)=0; in=in+1; %U 1
	nreads(in)=0; in=in+1; %V 2
	nreads(in)=0; in=in+1; %W 3
	nreads(in)=0; in=in+1; %TH 4
	nreads(in)=1; in=in+1; %Q01 5
	nreads(in)=1; in=in+1; %Q02 6
	nreads(in)=1; in=in+1; %Q03 7
	nreads(in)=1; in=in+1; %Q04 8
	nreads(in)=1; in=in+1; %Q05 9
	nreads(in)=1; in=in+1; %Q06 10
	nreads(in)=0; in=in+1; %Q07 11
	nreads(in)=0; in=in+1; %Q08 12
	nreads(in)=0; in=in+1; %Q09 13
	nreads(in)=0; in=in+1; %Q10 14 %tracer
	nreads(in)=1; in=in+1; %P 15 (pressure pert/rho)
    nreads(in)=1; in=in+1; %P 15 (pressure pert/rho)
    % N.B. might have to adjust these depending on no. of Q fields
end

%bezier=1;
%newton=0;

rezip=0; %flag say whether to rezip or not

zip=0;

lis=dir(direcDan(idir).dir);
for izip=1:length(lis)
    if strcmp(lis(izip).name,fname)==1 %if unzipped file is present then use that instead
        zip=0;
        break
    end
    if strcmp(lis(izip).name,[fname '.gz'])==1 %if zipped file present prepare to use that
        zip=1;
    end
end

if zip==1
    'unzipping file ...'
    eval(['!c:/cygwin/bin/gzip -d ''' FileName '.gz''']);
    'finished unzipping'
end

if rezip==0
    zip=0;  
end


if machine(idir)==2 %nlbas
    fid=fopen(FileName,'rb','s'); %Bezier/HPCx = 2 ('s') %give an array in machine to determine machine type - should match how are numbered in direcDan
    linux=1;
elseif machine(idir)==3
    fid=fopen(FileName,'rb','a'); %newton/Horace = 3
    linux=1;
elseif machine(idir)==4
    fid=fopen(FileName,'rb','l'); %Sami's Linux w/ recl=512
    linux=1;
elseif machine(idir)==5
    fid=fopen(FileName,'rb','l'); %Sami w/ Linux recl=128 (hopefully this will work!)
    linux=1;
elseif machine(idir)==6
    fid=fopen(FileName,'rb','l'); %Paul's Linux
    linux=0;    
elseif machine(idir)==7
    fid=fopen(FileName,'rb','s'); %HPCx 1024
    linux=2;    
elseif machine(idir)==8
    fid=fopen(FileName,'rb','a'); %
    linux=2;       
else
    fid=fopen(FileName,'rb'); %dos machine = 1
    linux=0;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just skip the header
%fseek(fid,18.*512.*4,'bof');

if linux==0|linux==2
    if linux==1
        ha=1024;
	else
        ha=2048;
	end
 
    
	hsize=512;
    

    if machine(idir)==7
        hsize=1024;
        ha=1024;
    end    
    
	ind=0;
	HEADER=((fread(fid,hsize,'2048*char')))';
	
	   
    
	
	
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


else %if linux==0

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

%if newton==1
% fseek(fid,1024+128*4,'cof');
%end

part=128;  %128
%APARMS=fread(fid,1024,'1024*float=>double');
%APARMS=fread(fid,1024,'float=>double');
[APARMS part]=Read1D_vrad2(1024,fid,linux,part);

% Now define the offsets
OFFSET(1)=(45.+APARMS(45)); %=ISTART1
OFFSET(2)=OFFSET(1) +double(APARMS(17));
OFFSET(3)=OFFSET(1) +double(2.*APARMS(17));
OFFSET(4)=OFFSET(1) +double(3.*APARMS(17));
OFFSET(5)=OFFSET(4) +double(49.);   %2nd ISTART2=ISTART+49+3*NQP
OFFSET(6)=OFFSET(5) +double(APARMS(17));
IHS=OFFSET(6) +double(15.);
NPARMS=IHS+67;

if ijustheader==1; %if just wanted to read in header & APARMS
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if APARMS(4)>1 & APARMS(1)==1
    NJ=APARMS(2)*APARMS(4); %parallel 2d case where no. of J points =JJP*NPES
else
    NJ=APARMS(2); %normal case where =JJP
end
ipep=ceil(APARMS(1)/APARMS(4)); %no. slices per processor

% if newton==1
% fseek(fid,1024*3*4,'cof');
% end


% if linux==1
% 	fseek(fid,4*1024,'cof');
% end

clear Grid

%fseek(fid,3*4*1024,'cof');

% In this version the blocks are blocks of 128
[Grid.X1 part]=Read1D_vrad2(APARMS(1)+2,fid,linux,part);
[Grid.X2 part]=Read1D_vrad2(APARMS(1)+2,fid,linux,part);


if APARMS(2)>1500  %needed to add for JJP=1536 run. not sure where threshold lies
%fseek(fid,1024*4,'cof');
end

[Grid.Y1 part]=Read1D_vrad2(NJ+2,fid,linux,part);
[Grid.Y2 part]=Read1D_vrad2(NJ+2,fid,linux,part);
[Grid.Z part]=Read1D_vrad2(APARMS(3),fid,linux,part);
[Grid.ZN part]=Read1D_vrad2(APARMS(3),fid,linux,part);



[Grid.RHO part]=Read1D_vrad2(APARMS(3),fid,linux,part);


[Grid.RHON part]=Read1D_vrad2(APARMS(3),fid,linux,part);


% Choice of IDGx dependent
if(APARMS(OFFSET(5)+1) > 0)
    [Grid.UBAR part]=Read1D_vrad2(APARMS(3),fid,linux,part);
end
if(APARMS(OFFSET(5)+2) > 0)
    [Grid.VBAR part]=Read1D_vrad2(APARMS(3),fid,linux,part);
end
if(APARMS(OFFSET(5)+5) > 0)
    [Grid.OLTHBAR part]=Read1D_vrad2(APARMS(3),fid,linux,part);
    [Grid.THINIT part]=Read1D_vrad2(APARMS(3),fid,linux,part);
    [Grid.THREF part]=Read1D_vrad2(APARMS(3),fid,linux,part);
end
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) > 0)
        [Grid.OLQBAR(1:APARMS(3),IQ) part]=Read1D_vrad2(APARMS(3),fid,linux,part);
    end
end
if(APARMS(OFFSET(5)+6) > 0)
    [Grid.PREFN part]=Read1D_vrad2(APARMS(3),fid,linux,part);
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

     [SENflux part]=Read1D_vrad2(NJ,fid,linux,part);
     SENflux=SENflux.*Grid.RHO(1).*1005; %convert to W/m2
     [LATflux part]=Read1D_vrad2(NJ,fid,linux,part)
     LATflux=LATflux.*Grid.RHO(1).*2.501e6;
 end
 
 if APARMS(NPRAB+2)==1 % read in surface varying roughness lengths

     [Z0var part]=Read1D_vrad2(NJ,fid,linux,part);
     [Z0THvar part]=Read1D_vrad2(NJ,fid,linux,part);
 end



% Now XXX
if(APARMS(35) >= 1 & APARMS(36) ~= 0)
    % XXX_A
    for IL=1:APARMS(35)
        [Grid.XXX_A(1:APARMS(3),IL) part]=Read1D_vrad2(APARMS(3),fid,linux,part);
    end
    
    % XXX_W
    for IL=1:APARMS(35)
        [Grid.XXX_W(1:APARMS(3),IL) part]=Read1D_vrad2(APARMS(3),fid,linux,part);
    end
        
    % XXX_TH
    for IL=1:APARMS(35)
        [Grid.XXX_TH(1:APARMS(3),IL) part]=Read1D_vrad2(APARMS(3),fid,linux,part);
    end
        
    % XXX_QYY
    for IQ=1:APARMS(17)
        for IL=1:APARMS(35)
            [Grid.XXX_QYY(1:APARMS(3),IQ,IL) part]=Read1D_vrad2(APARMS(3),fid,linux,part);
        end
    end
end
% Import the time-averaged diagnostics
if(APARMS(OFFSET(6)+8) ==1)
    [TimeAv.RNDGS part]=Read1D_vrad2(APARMS(OFFSET(6)+9),fid,linux,part);
    [TimeAv.DGAV part]=Read2D_vrad2(APARMS(3),APARMS(OFFSET(6)+9),fid,linux,part);
    TimeAv.DGAV=TimeAv.DGAV';
    %TimeAv.DGAV=Read1D(APARMS(OFFSET(6)+9).*APARMS(3),fid);
end


if(APARMS(NPRAB+11) >0)
    for k=1:APARMS((NPRAB+11))
%        [TimeAv.surav(:,:,k) part]=Read2D_vrad2(NJ,APARMS(1),fid,linux,part);
%        [TimeAv.surav(:,:,k) part]=Read1D_vrad2(APARMS(1),fid,linux,part);
        [TimeAv.surav(:,:,k) part]=Read2D_vrad2(ipep,APARMS(2),fid,linux,part);
    end
end



if(APARMS(OFFSET(6)+14) ==1 & APARMS(OFFSET(6)+13) >=1)
    %Grid.DGSPAV=Read1D2(APARMS(OFFSET(6)+13),floor(JJTOTP./2+1),fid,linux,part);
    [Grid.DGSPAV part]=Read1D_vrad2(APARMS(OFFSET(6)+13).*(APARMS(44)./2+1),fid,linux,part);
end
% Time Series
if(APARMS(28) == 1)
   for NSER=1:APARMS(34)
        [X,part]=Read1D_vrad2(APARMS(OFFSET(6)+10),fid,linux,part);
        SER(:,NSER)=X;
    end
end

if (APARMS(37)==0) %if no radiation
    SER(:,38:APARMS(34)+7)=SER(:,31:APARMS(34));  %move timeseries down so that process rates are in the places listed on sheet
    SER(:,31:37)=0;
end

nqp=APARMS(17);
diff=nqp-9;

nprates=70;

%maxq01 goes in slot 38
%hq01 therefore goes in 38+nqp and then PUD01 in 38+2*nqp

%NEXT=nts_pudxx+nmetp passed to DansTimSer for the next slot in SER = 38+2*nqp + 4
%first values (prates) entered in NEXT+iq-1 where iq=1:nprates
%radmax therefore goes in NEXT+nprates
%nts_pudxx is the first slot for the PUDxx values (PUD01)

%NEXT=38+2*nqp+4
%then radmax is in 38+2*nqp+4+nprates

endno=2*imaxrad+itrtop+2*imaxpr+2*icape+2*iaero;
%col 108=maxrad 109=maxrad height 110=max tracer height
%col 111=max precipitation rate 112=J location of max pp rate
%col 113=average cape, 114=av CIN
%115=max supersat (ev/es) (Nenes scheme), 116=height of max

imphys=38+2*nqp + 4;
%%%%%%%%%%%% gets process rate names and puts them in correct order into pname(i).p
pnames; %%%% PGDEP, etc




iadjust=0;
if iadjust==1
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
    [TwoD.U part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
end
if(APARMS(OFFSET(5)+2) == 2)
    [TwoD.V part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
end
if(APARMS(OFFSET(5)+3) == 2)
    [TwoD.W part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
end
if(APARMS(OFFSET(5)+5) == 2)
    [TwoD.TH1 part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
end


% 2-D Q variables
for IQ=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) == 2)
      %  [TwoD.Q(:,:,IQ) part]=Read2D_vradNOSKIP(NJ+2,APARMS(3),fid,linux,part);
                [TwoD.Q(:,:,IQ) part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);

    end
end









% RADIATION HERE ***********************
if(APARMS(37)>=1)
    [TwoD.sthrad part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
    if(APARMS(38)==1)
        [TwoD.twdrad part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
    end
    if(APARMS(39)==1)
        [TwoD.swdrad part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
    end
else
    TwoD.sthrad=0;
    TwoD.twdrad=0;
    TwoD.swdrad=0;
end

% more code here
   % if(APARMS(37)>=2)
   %     if(APARMS(38)==1)

%try      %%%%%%%%%%%%%%%%%
    
if(APARMS(OFFSET(5)+6) == 2)
    [TwoD.P part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
    RHOref=repmat(Grid.RHON,[1 length(Grid.Y1)]);
    Pref=repmat(Grid.PREFN,[1 length(Grid.Y1)]);
    TwoD.PP=TwoD.P.*RHOref + Pref;
else
    TwoD.P=0;
    TwoD.PP=0;
end
if(APARMS(OFFSET(5)+7) == 2)
    [TwoD.PV part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
else
    TwoD.PV=0;
end
if(APARMS(40) == 2)
    [TwoD.TH2 part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
else
    TwoD.TH2=0;
end
if(APARMS(41) == 2)
    [TwoD.THV part]=Read2D_vrad2(NJ+2,APARMS(3),fid,linux,part);
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
    %TEMP=Read1D_vrad2(APARMS(3),fid);
    
if(APARMS(NPRAB+13)==0) %igaloffp=0 - Galilean transform
    [surfdiag.precgal part]=Read2D_vrad2(APARMS(NPRAB+16),APARMS(NPRAB+15),fid,linux,part); 
    [surfdiag.pmaxgal part]=Read2D_vrad2(APARMS(NPRAB+16),APARMS(NPRAB+15),fid,linux,part); 
    [surfdiag.instantgal part]=Read2D_vrad2(APARMS(NPRAB+16),APARMS(NPRAB+15),fid,linux,part); 
end
    

if(APARMS(NPRAB+14)==1)
    if APARMS(1)>1
        if APARMS(3)>1
            nis=APARMS(1)+2;
            njs=APARMS(2)+2;
        else
            %nis=APARMS(1);
        end
    else
        njs=APARMS(2);
        nis=1;
    end
    
 %   njs=APARMS(2); %may need to change this to NJ if the 2-D parallel case is changed so that it outputs correctly (JJP*NPES)
                    %instead of just JJP arrays.
    
    
    for imetp=1:APARMS(45)
	    if(APARMS(45+imetp)==2 | APARMS(45+imetp)==3)
            [surfdiag.prec(1:nis,1:njs,imetp) part]=Read2D_vrad2(njs,nis,fid,linux,part);    %total acc precip over ALL of runtime
	    end
    end
    [surfdiag.pmax part]=Read2D_vrad2(njs,nis,fid,linux,part); %max prate over diag time
    [surfdiag.class part]=Read2D_vrad2(njs,nis,fid,linux,part); 
    [surfdiag.citop part]=Read2D_vrad2(njs,nis,fid,linux,part); 
    [surfdiag.lwp part]=Read2D_vrad2(njs,nis,fid,linux,part);
    [surfdiag.iwp part]=Read2D_vrad2(njs,nis,fid,linux,part);
    [surfdiag.instant part]=Read2D_vrad2(njs,nis,fid,linux,part);
end


xlen=NJ+2;

iread=1; %counter that gets incremented with every call to Read3D2

if(APARMS(OFFSET(5)+1) == 3)
%      for i=1:xlen^2
%          fseek(fid,8,'cof');
%          fseek(fid,4,'cof');
%     end
%     
%     for j=1:xlen
%             fseek(fid,8,'cof');
%             fseek(fid,8,'cof');
%             fseek(fid,8,'cof');
%             fseek(fid,8,'cof');
%             fseek(fid,8,'cof');
%             fseek(fid,8,'cof');
%     end#
    
    [ThreeD.U part iread]=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);

end



if(APARMS(OFFSET(5)+2) == 3)
   [ThreeD.V part iread]=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);
end

if(APARMS(OFFSET(5)+3) == 3)
    [ThreeD.W part iread]=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);
end

if(APARMS(OFFSET(5)+5) == 3)
    [ThreeD.TH1 part iread]=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);
end


% 3-D Q variables
if idontclose==1
    return
end

istore2=0;
for istore=1:APARMS(17)
    if(APARMS(OFFSET(5)+7+IQ) == 3)
        if nreads(iread)==1
        	istore2=istore2+1;    
            [ThreeD.Q(:,:,:,istore2) part iread]=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);
        else 
            [temp part iread]=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);
        end                
%            ThreeD.Q(:,:,:)=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);

    end
end    

if(APARMS(OFFSET(5)+6) == 3)
    [ThreeD.P part iread]=Read3D2(APARMS(1)+2,NJ+2,APARMS(3),fid,linux,part,iread,nreads);
%    RHOref=repmat(Grid.RHON,[1 length(Grid.Y1)]);
%    Pref=repmat(Grid.PREFN,[1 length(Grid.Y1)]);
%    TwoD.PP=TwoD.P.*RHOref + Pref;
elseif (APARMS(1) > 1) %IF 3D
    ThreeD.P=0;
    ThreeD.PP=0;
end


% catch
%     fprintf(1,'WARNING PROBLEM WITH FILE jj=%d',jj); %because in ccn960 newt case is some problem with file 17
% end

%c=ftell(fid)

if idontclose==0

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


if zip==1
    're-zipping file ...'
    command=['!c:/cygwin/bin/gzip ''' FileName ''''];
    eval(command);
end

end