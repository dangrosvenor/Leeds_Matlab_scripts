function [m,dz]=TotMassBudgetProf(Grid,diag,time,iz,iz2)

%Grid.Z=[-75;Grid.Z];

%diag(iz+1:iz2,:)=diag(iz+1:iz2,:);
%diag(iz,:)=0;

%[dfdz dz]=ddz(Grid.Z(iz:iz2),diag(iz:iz2,1:end)); %finds d/dz of flux
%dt=repmat(time(2:end)-time(1:end-1),[iz2-iz 1])*3600; %dt in secs

dz=Grid.Z(2:end)-Grid.Z(1:end-1);

rho=repmat(Grid.RHO(iz-1:iz2-1),[1 size(diag,2)]);

%mtim=dfdz.*dz.*rho.*dt;

%mtim=dfdz.*dz.*rho.*300;


%mtim=dfdz.*dt;

%m=sum(mtim,2);

tzc1=repmat(Grid.RHO(iz-1:iz2-2)./Grid.RHON(iz-1:iz2-2)./dz(iz:iz2-1) , [1 size(diag,2)] );
tzc2=repmat(Grid.RHO(iz:iz2-1)./Grid.RHON(iz:iz2-1)./dz(iz:iz2-1) , [1 size(diag,2)] );


%[df]=sum(diag(iz+1:iz2,:).*rho(2:end,:) - diag(iz:iz2-1,:).*rho(1:end-1,:) ,2)./dz;

[df]=diag(iz+1:iz2,:).*tzc2 - diag(iz:iz2-1,:).*tzc1;

m= df*300;