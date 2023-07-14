[iz,iz2]=findheight(GridDan(1).Z+620,17e3,20e3);
[t1]=findheight(GridDan(1).t+3,23.5);
t2=72; %final dump

microicerate=f*sum(icediagsALL(idir).i(iz:iz2,t1:t2,31:33),3)/npes; %ice mixing ratio source rate
profile=sum(microicerate,2);

rho=GridDan(idir).RHON(iz:iz2);

sumq=cumsum(flipud(profile).*rho);