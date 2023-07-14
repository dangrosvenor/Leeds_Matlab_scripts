idir=1;

H=16.15e3;


gamma=9.8e-3; %dry adiabatic lapse rate
g=9.8; % gravity

T=TempLes(GridDan(idir));
z=GridDan(idir).Z+620;

ih=findheight(z,H);

dtdz=diff(T)./diff(z);


N=sqrt( g./T(2:end) .* (dtdz + gamma) );

Nz=N(ih)


rho=GridDan(idir).RHON;
drhodz=diff(rho)./diff(z);

N2=sqrt(-g./rho(2:end) .* drhodz);
N2(ih)


dlnthdz=diff(log(GridDan(idir).THREF)) ./ diff(z);

N3= sqrt (g * dlnthdz);
N3(ih)

