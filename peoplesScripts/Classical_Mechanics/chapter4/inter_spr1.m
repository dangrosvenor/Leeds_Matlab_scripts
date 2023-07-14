%inter_spr1.m
%plots the coordinate solutions for the single mode coupled spring-mass 
%system without walls
clear;
m1=1;                                %masses
m2=2;
k0=0.5;                              %spring constant
x10=1;                               %initial positions
x20=-1;
v10=0.02;                            %initial speeds
v20=0.04;
mu=m1*m2/(m1+m2);                    %reduced mass
xcm0=(m1*x10+m2*x20)/(m1+m2);        %initial center of mass
vcm=(m1*v10+m2*v20)/(m1+m2);         %center of mass speed
xr0=(x20-x10);                       %relative coordinate
vr0=(v20-v10);                       %relative speed
om=sqrt(k0/mu);                       %frequency
tau=2*pi/om;                         %period
A=vr0/om;                            %amplitudes
B=xr0;
tmax=2;
str=cat(2,'m_1=',num2str(m1),',m_2=',num2str(m2),',k_0=',num2str(k0),...
',x_{10}=',num2str(x10),',x_{20}=',num2str(x20),',v_{10}=',num2str(v10),...
',v_{20}=',num2str(v20),',tmax=',num2str(tmax),'\tau');
tmax=tmax*tau;
t=[0:tau/50:tmax];                   %plotting time range
xr=A*sin(om*t)+B*cos(om*t);          %solution
xcm=xcm0+vcm*t;                      %cm position vs time
x1=xcm-m2*xr/(m1+m2);                %mass positions versus time
x2=xcm+m1*xr/(m1+m2);
plot(t,xcm,'k-.',t,x1,'b-',t,x2,'r--');
axis([0 tmax -1.5  1.5]);
text(0.25,1.35,str,'FontSize',11,'Color','black');
title('Single Mode Spring-Mass System Without Walls','FontSize',14)
h=legend('xcm','x1','x2',4); set(h,'FontSize',14)
ylabel('Position','FontSize',14);
xlabel('Time','FontSize',14);
