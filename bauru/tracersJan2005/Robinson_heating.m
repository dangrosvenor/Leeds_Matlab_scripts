A0=3e-8;
b=0.1;
Lx=600;
x=[0:0.1:Lx]/Lx;
T=273+30;

p=A0*exp(-(x-0.5).^2/b)*T;

%figure;
%plot(x*Lx,p);

dt=3;
t=[0:dt:3600];
dTdt=T*A0*t/dt;
dT=sum(dTdt)*dt;

figure;
plot(t,dTdt)