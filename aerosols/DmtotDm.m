function DmtotDm=lognormal(m,N,sig,Dg,rhoS)
%returns dMtot/dM for lognormal distribution of total number N standard deviation sig and number mode diameter Dg
%for calculation of total mass of mass bin
D=(m*3/4/pi/rhoS*1000).^(1/3)*2/100; 

DmtotDm = 1/3 * N/sqrt(2*pi)/log(sig) * exp( - (log(D/Dg)).^2/2/(log(sig)^2) );