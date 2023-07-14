clear wav wav2 dzz sw


for idir=1:3

dq=1000*(dq_tot(idir).d(:,:,1)...
         - nn(idir).n(:,:)/502 .* (3.67 - 3.8 ) );
     
     
z=GridDan(idir).Z+620;
dzz=(z(2:end)-z(1:end-1))/1000;
dzz(end+1)=dzz(end);
     
[z0 zend]=findheight(z,15e3,20e3);

dz=0;
ndz=1;
while (dz<1e3 & ndz<length(z)+1)
    dz=dz+z(z0+ndz+1)-z(z0+ndz);
    ndz=ndz+1;
end

for it=1:size(dq,2)
for idc=z0:zend-ndz+1
    winav(idc-z0+1)=mean(dq(idc:idc+ndz-1,it));
end

wav(idir).w(it)=max(winav);
wav2(idir).w(it)=sum(dq(:,it).*GridDan(idir).RHON.*dzz);
end


dw=wav2(idir).w(2:end)-wav2(idir).w(1:end-1);
idw=find(dw<0);
dw(idw)=0;
sw(idir).s=cumsum(dw);
end

