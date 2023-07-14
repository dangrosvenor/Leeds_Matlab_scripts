function [m,dz]=TotMassBudgetProfRate(Grid,diag,time,iz,iz2)
%takes a dq/dt diag matrix and calculates time integrated profile (*rho)
dz=diffrep(Grid.Z(iz:iz2),diag(iz:iz2,1:end)); %finds d/dz of flux
%dt=repmat(time(2:end)-time(1:end-1),[iz2-iz 1])*3600; %dt in secs

rho=repmat(Grid.RHO(iz:iz2-1),[1 length(time)]);


%dt=300 seconds;

m=sum(diag(iz+1:iz2,1:end).*dz.*rho,2)*300;



