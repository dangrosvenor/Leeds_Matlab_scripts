function wi=LEWW(T,P)
T0=273.15;
c1=3.8;
c2=-17.2693882;
c3=35.86;
c4=6.109;

P=P./100;
wi=c1./(P.*exp(c2.*((T-T0)./(T-c3)))-c4);