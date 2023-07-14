clear m

% [dfdz(1).d dz]=ddz(Grid.Z,icediag3(1).i(:,1:end-1,10)); %finds d/dz of flux
% dt=time(2:end)-time(1:end-1);
% 
% rho=repmat(Grid.RHO(2:end),[1 length(time)-1]);
% 
% mtim=sum(dfdz(1).d.*dz.*rho,1);
% 
% m(1).m=sum(mtim.*dt);


if ~exist('icediag3')
load('g:runs/dmi1715_5ppmv/results/diags/icediags3');
end
if ~exist('icediag4_2') %ones that are not divided by area
load('g:runs/dmi1715_5ppmv/results/diags/icediags4_2');
end
if exist('Grid')~=1
	FileName='g:runs/dmi1715_5ppmv/RUN0001.DG0001';
	DIAG_3d;
end
if exist('grid')~=1
load('g:runs/dmi1715_5ppmv/results/diags/grid');
end



[iz,iz2]=findheight(Grid.Z,0e3,30e3);
iz2=length(Grid.Z);
iz=2;

%see note3 for variables in icediag3
i=1;
m(1).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,10)+icediag3(i).i(:,:,16),time,iz,iz2);
%m(2).m=sum(sum(icediag3(i).i(iz+1:iz2,2:end,14).*rho(iz+1:iz2,:).*dz(iz+1:iz2,:),1 ).*dt);
m(3).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,3),time,iz,iz2); %FQ06
m(4).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,6),time,iz,iz2); %FQ04
m(5).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,9),time,iz,iz2); %FQ05
m(6).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,1),time,iz,iz2); %WQ06
m(7).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,4),time,iz,iz2); %WQ04
m(8).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,7),time,iz,iz2); %WQ05



m(9).m=m(3).m+m(4).m+m(5).m; %sum of fqs
m(10).m=m(6).m+m(7).m+m(8).m; %sum of wqs for ice


ovall=-m(1).m-m(10).m+m(9).m;

m(11).m=TotMassBudgetProfRate(Grid,icediag3(i).i(:,:,14),time,iz,iz2); %from ALL_DQ01

m(12).m=TotMassBudgetProfRate(Grid,icediag3(i).i(:,:,12),time,iz,iz2); %from ALL_DQ06
m(13).m=TotMassBudgetProfRate(Grid,icediag3(i).i(:,:,13),time,iz,iz2); %from ALL_DQ04
m(14).m=TotMassBudgetProfRate(Grid,icediag3(i).i(:,:,15),time,iz,iz2); %from ALL_DQ05

m(15).m=m(12).m+m(13).m+m(14).m;

m(16).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,10),time,iz,iz2); %ALL_WQ01
m(17).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,16),time,iz,iz2); %ALL_WQSG01
m(18).m=TotMassBudgetProf_mass(Grid,icediag3(i).i(:,:,17),time,iz,iz2); %ALL_WQSG04

for i=1:size(icediag4(1).i,3)
    m(18+i).m=TotMassBudgetProf_mass(Grid,icediag4(1).i(:,:,i),time,iz,iz2); 
end

%m(9).m=m(3).m+m(4).m+m(5).m; %sum of fqs
%m(10).m=m(6).m+m(7).m+m(8).m + m(18).m + m(18+11).m + m(18+12).m + m(18+15).m + m(18+16).m; %sum of wqs for ice

totw=sum(icediag3(1).i(:,:,[1 4 7 10 16]),3);

m(37).m=TotMassBudgetProf_mass(Grid,totw,time,iz,iz2);

toticesub=sum(icediag4(1).i(:,:,[5:16]),3); %total water flux with sub grids
m(38).m=TotMassBudgetProf_mass(Grid,toticesub,time,iz,iz2);

