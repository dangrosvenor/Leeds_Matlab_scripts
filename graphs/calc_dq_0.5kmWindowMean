dq=1000*(dq_tot(idir).d(izmin:izmax,dumprange,1)...
         - nn(idir).n(izmin:izmax,dumprange)/502 .* (3.67 - 3.8 ) );
     
     
z=GridDan(idir).Z;
     
[z0 zend]=findheight(z,15e3,20e3);

dz=0;
ndz=1;
while (dz<0.5e3 & ndz<length(z)+1)
    dz=dz+z(z0+ndz+1)-z(z0+ndz);
    ndz=ndz+1;
end

for it=1:size(dq,2)
for idc=z0:zend-ndz+1
    winav(idc-z0+1)=mean(dq(idc:idc+ndz-1,it));
end

wav(it)=max(winav);
end