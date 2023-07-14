%dis=[0.01:0.01:0.5];
dis=0.5;
N=1000e6;
ql=[0:0.1:10].*1e-3; %[g/kg]*1e-3 = kg/kg
rhoair=1.083;


qt=pi*20e-6^3*1000*N/6/rhoair;



alx=dis.^-2 - 1;

b6=(gamma(6+alx+1)./gamma(1+alx)).^(1/6) .* (gamma(3+alx+1)./gamma(1+alx)).^(-1/3);

%figure;
%plot(dis,b6.^6);


liu=1.9e17*b6^6*(3*rhoair/4/pi/1000)^2/N *(ql).^3;
kess=10e-4*(ql-qt);

ii=find(ql<qt);

liu(ii)=0;
kess(ii)=0;

figure;
plot(ql*1000,liu,'k-');
hold on;
plot(ql*1000,kess,'k.-');

