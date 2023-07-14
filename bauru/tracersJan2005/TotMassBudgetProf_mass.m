function [m,dz]=TotMassBudgetProf(Grid,diag,time,iz,iz2)

%Grid.Z=[-75;Grid.Z];

%diag(iz+1:iz2,:)=diag(iz+1:iz2,:);
%diag(iz,:)=0;

%[dfdz dz]=ddz(Grid.Z(iz:iz2),diag(iz:iz2,1:end)); %finds d/dz of flux
%dt=repmat(time(2:end)-time(1:end-1),[iz2-iz 1])*3600; %dt in secs



rho=repmat(Grid.RHO(iz-1:iz2-1),[1 length(time)]);

%mtim=dfdz.*dz.*rho.*dt;

%mtim=dfdz.*dz.*rho.*300;


%mtim=dfdz.*dt;

%m=sum(mtim,2);


[df]=diag(iz+1:iz2,:).*rho(2:end,:) - diag(iz:iz2-1,:).*rho(1:end-1,:);

m= sum(df,2)*300;