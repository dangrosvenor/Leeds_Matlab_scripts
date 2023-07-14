%foft.m
clear;
m=1.0;                             %mass
f0=1.0;                            %force amplitude
w=3.0;                             %force angular frequency
x0=0.0;                            %initial position
v0=0.05;                           %initial velocity
t=[0:0.1:10];                      %time array
a=f0*cos(w.*t)/m;                  %acceleration array
v=v0+f0*sin(w.*t)/m/w;             %velocity array
x=x0+v0.*t-f0*(cos(w.*t)-1)/m/w/w; %displacement array
plot(t,x,'k-',t,v,'b:',t,a,'r-.');
title('x(t),v(t),a(t) due to F=F_0cos(\omegat)','FontSize',14)
ylabel('x(t),v(y),a(t)','FontSize',14);
xlabel('time(sec)','FontSize',14);
h=legend('x(t)','v(t)','a(t)'); set(h,'FontSize',14)