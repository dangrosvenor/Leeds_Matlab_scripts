%dis=[0.01:0.01:0.5];
dis=0.6;
N=[240 480 720];
N=N*1e6;
rhoair=1.083;
ql=[0:0.1:4].*1e-3/rhoair; %[g/kg]*1e-3 = kg/kg
figure;

for ic=1:size(N,2)

qt=pi*20e-6^3*1000*N(ic)/6/rhoair;



alx=dis.^-2 - 1;

b6=(gamma(6+alx+1)./gamma(1+alx)).^(1/6) .* (gamma(3+alx+1)./gamma(1+alx)).^(-1/3);

%figure;
%plot(dis,b6.^6);


liu=1.9e17*b6^6*(3*rhoair/4/pi/1000)^2/N(ic) *(ql).^3;
kess=10e-4*(ql-qt);

ii=find(ql<qt);

liu(ii)=0;
kess(ii)=0;

subplot(size(N,2),1,ic);

plot(ql*1000*rhoair,liu,'b-x');
hold on;
plot(ql*1000*rhoair,kess,'k-o');
titl=strcat('dis=',Num2str(dis),',N=',Num2str(N(ic)/1e6),' /cm^3');
title(titl);

end

