function [m,dz]=TotMassBudget(Grid,diag,time,iz,iz2)

[dfdz dz]=ddz(Grid.Z(iz:iz2),diag(iz:iz2,1:end-1)); %finds d/dz of flux
dt=time(2:end)-time(1:end-1);

rho=repmat(Grid.RHO(iz+1:iz2),[1 length(time)-1]);
function [m,dz]=TotMassBudget(Grid,diag,time,iz,iz2)

[dfdz dz]=ddz(Grid.Z(iz:iz2),diag(iz:iz2,1:end-1)); %finds d/dz of flux
dt=time(2:end)-time(1:end-1);

rho=repmat(Grid.RHO(iz+1:iz2),[1 length(time)-1]);

mtim=sum(dfdz.*dz.*rho,1);

m=sum(mtim.*dt);function [m,dz]=TotMassBudget(Grid,diag,time,iz,iz2)

[dfdz dz]=ddz(Grid.Z(iz:iz2),diag(iz:iz2,1:end-1)); %finds d/dz of flux
dt=time(2:end)-time(1:end-1);

rho=repmat(Grid.RHO(iz+1:iz2),[1 length(time)-1]);

mtim=sum(dfdz.*dz.*rho,1);

m=sum(mtim.*dt);
mtim=sum(dfdz.*dz.*rho,1);

m=sum(mtim.*dt);





