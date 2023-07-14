%radar calculation checking

iz=100;
ix=50;

[lam,nx0,ns,ms,D]=gamlemRow ( TwoD.Q(iz,ix,9), TwoD.Q(iz,ix,4) , GridDan(1).RHON(iz),'s');
dD=diff(D);
radar=10*log10( 0.224*sum(dD.*ns(2:end).*(D(2:end)*1e3).^6) ); %convert D to mm as Z units in mm6/m3 ns=dN/dD /m3 /m so *D


a=find(D<1e-2);
a=a(end);
%dD=diff(D(1:a));
%radar=10*log10( 0.224*sum(dD.*ns(2:a).*(D(2:a)*1e3).^6) ); %convert D to mm as Z units in mm6/m3 ns=dN/dD /m3 /m so *D


