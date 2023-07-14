dz=z(2:end)-z(1:end-1);
d=size(icediagALu(1).i(2:end,:,5),2);
dz=repmat(dz,[1 d]);
% wq=[zeros([1 64]) ; ( icediagALu(1).i(2:end,:,5) - icediagALu(1).i(1:end-1,:,5) )./dz];
wq=[zeros([1 64]) ; ( icediagALu(1).i(2:end,:,5) )./dz];

sumi=[3 4 5 6 7 9 11 12 13];

		dgs{1}='PIMLT';
		dgs{2}='PSAUT';
		dgs{3}='PSACI';
		dgs{4}='PRACI_S';
		dgs{5}='PGACI';
		dgs{6}='PRACI_G';
		dgs{7}='PIHAL';
		dgs{8}='PIPRM';
		dgs{9}='PICNT';
		dgs{10}='PIDEP';
		dgs{11}='PIACW';
		dgs{12}='PIFRW';
        dgs{13}='RSAUT';
        dgs{14}='RIACI';
        dgs{15}='PISUB';
        
		dgs{16}='PSDEP';
		dgs{17}='PIACR_S';
		dgs{18}='PSACW';
		dgs{19}='PSSUB';
		dgs{20}='PGACS';
		dgs{21}='PRACS';
		dgs{22}='PGAUT';
		dgs{23}='PSMLT';
		dgs{24}='RSBRK';
		dgs{25}='RIACR_S';
		dgs{26}='RGACS';
        dgs{27}='RSACR';
        dgs{28}='RSACS';
        

%for calc of dqNI        
MI0=1e-15;
imi0=[7 8 9 12]; %indices for processes to be divided by pristine ice mass MI0- PIHAL,PICNT,PIPRM,PIFRW to get the number source
imnq=[3 4 5 6]; %to be multiplied by n/q to get the number conc
inorm=[13 14]; %no factor needed


        
dqtot=icediag2(1).i(:,:,4) - ( sum(icediag(1).i(:,:,imi0),3)/MI0 - sum(icediag(1).i(:,:,inorm),3) );
%difference between dq from average dq from diags and that from process rates to get approx value of ones that need
%to be multiplied by n/q.
%should also take into account the effect of the noumber conc being limited by RNc - can't really do both at the same time.

dmnq=sum(icediag(1).i(:,:,imnq),3);
ndivqav= - dqtot./dmnq; %av value of n/q to match dq from diag

ii=find(ndivqav==-Inf);
ndivqav(ii)=0;
jj=find(ndivqav<0);
ndivqav(jj)=0; %check how many points lead to negative n/q

ndivqav(dqtot<(0.001))=0;
ndivqav(dmnq<(0.001))=0; %set to zero in situations where either dq or dmnq are near zero to avoid spurious values of n/q

ndivqav(ndivqav==Inf)=0; %put Inf values as zeros

sumLiq= sum(icediag(1).i(:,:,imnq),3).*ndivqav + sum(icediag(1).i(:,:,[7]),3)/MI0 + icediag(1).i(:,:,13);

sumLiqNonMPC=wq.*sumLiq;

ii=find(ndivqav==-Inf);
ndivqav(ii)=0;
jj=find(ndivqav<0);
ndivqav(jj)=0; %check how many points lead to negative n/q


%for calc of qs/ns       

imnq=[19 23 22]; %to be multiplied by n/q
inorm=[25 13 24]; %no factor needed
inormneg=[26 27 28]; %neagtive values

%ALL_DQ09 = icediag2(1).i(:,:,10)        
dqtot=icediag2(1).i(:,:,10) - sum(icediag(1).i(:,:,inorm),3) + sum(icediag(1).i(:,:,inormneg),3);
%difference between dq from average dq from diags and that from process rates to get approx value of ones that need
%to be multiplied by n/q.

dqnq=sum(icediag(1).i(:,:,imnq),3);
ndivqavSnow= - dqtot./dqnq; %av value of n/q to match dq from diag

ndivqavSnow(dqtot<(0.001))=0; %set to zero when either dq or dqnq are close to zero to avoid spurious values of n/q
ndivqavSnow(dqnq<(0.001))=0;

ii=find(ndivqavSnow==-Inf);
ndivqavSnow(ii)=0;
jj=find(ndivqavSnow<0);
ndivqavSnow(jj)=0; %check how many points lead to negative n/q
ndivqavSnow(ndivqavSnow==Inf)=0; %put Inf values as zeros
