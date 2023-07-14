wq=icediag2(1).i(:,:,5);

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

%for calc of dqNI        
MI0=1e-15;
imi0=[7 8 9 12]; %indices for processes to be multiplied by MI0
imnq=[3 4 5 6]; %to be multiplied by n/q
inorm=[13 14]; %no factor needed
        
dqtot=icediag2(1).i(:,:,4) - ( sum(icediag(1).i(:,:,imi0),3)*MI0 - sum(icediag(1).i(:,:,inorm),3) );
%difference between dq from average dq from diags and that from process rates to get approx value of ones that need
%to be multiplied by n/q.

ndivqav=dqtot./(sum(icediag(1).i(:,:,imnq),3)); %av value of n/q to match dq from diag
ndivqav(ndivqav==Inf)=0;
sumLiq= sum(icediag(1).i(:,:,imnq),3).*ndivqav + sum(icediag(1).i(:,:,[7 9 12]),3)*MI0 - icediag(1).i(:,:,13);



%follows are the labels for the remaining ice and snow d/dt calcs stored in icediagsSnow(1).i

dgs{1}='PISUB';
dgs{2}='PSDEP';
dgs{3}='PIACR_S';
dgs{4}='PSACW';
dgs{5}='PSSUB';
dgs{6}='PGACS';
dgs{7}='PRACS';
dgs{8}='PGAUT';
dgs{9}='PSMLT';
dgs{10}='RSBRK';
dgs{11}='RIACR_S';
dgs{12}='RGACS';
dgs{13}='RSACR';
dgs{14}='RSACS';

















