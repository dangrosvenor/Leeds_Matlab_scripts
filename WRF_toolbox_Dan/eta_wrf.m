function [eta,zeta,p,psfc]=eta_wrf(nc)


psfc=nc{'PSFC'}(1,:);
psfc=psfc(1,1);
p=WRFUserARW(nc,'P',1,1,1);
p=p*100;
eta=([psfc p]-5000)/(psfc-5000); %5000 is the pressure top
zeta=[0 WRFUserARW(nc,'Z',1,1,1)];