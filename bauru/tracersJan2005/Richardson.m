%estimates the Richardson number at each 2-d grid point
%Ri=N^2/(dV/dz)^2 where N^2=g/th dth/dz (N=Brunt-Vas frequency, th=potemp)

clear diff

z=GridDan(1).Z+620;
dzz=diff(z);
dz=repmat(dzz,[1 size(TwoD.TH2,2)]);

th=TwoD.TH2; %potemp
dthdz=diff(th)./dz;
dVdz=diff(TwoD.V)./dz;
N2=9.81./th(1:end-1,:) .* dthdz;
Ri=N2./(dVdz).^2;

'done Richardson'

