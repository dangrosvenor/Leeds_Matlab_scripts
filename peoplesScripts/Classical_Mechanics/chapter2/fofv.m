%fofv.m
clear;
g=input(' Enter value of gravity g: ');
m=input(' Enter the object mass m: ');
C=input(' Enter the drag coefficient C: ');
if C< 1.e-3 C=1.e-3; end                     %prevent division by zero
y0=input(' Enter the initial height y0: ');  %Initial Conditions
v0=input(' Enter the initial velocity v0: ');%Initial Conditions
%if needed use next line to make NPTS an input
%NPTS=input(' Enter the number calculation steps desired NPTS: ');
NPTS=100;
% use next line to make TMAX a desired input
%TMAX=input(' Enter the run time TMAX: ');
%tz is the time for the frictionless free fall case to have y=0
tz=v0/g+sqrt((v0/g)^2+2*y0/g);
%define a function of several parameters using the inline method
f=inline('y0-m*(g*t+(m*g/C+v0)*(exp(-C*t/m)-1))/C','t','g','m','C','y0','v0');
%TMAX is the time for the case with friction to reach zero, use tz as guess
TMAX=fzero(f,tz,[],g,m,C,y0,v0);
t=[0:TMAX/NPTS:TMAX];%time array from zero to the time to reach ground
y=y0-m*(g*t+(m*g/C+v0)*(exp(-C*t/m)-1))/C;
v=(m*g/C+v0)*exp(-C*t/m)-m*g/C;
a=-g-C*v/m;
plot(t,y,'k-',t,v,'b:',t,a,'r-.');
title('Resistive Force Proportional to -v Example','FontSize',14)
ylabel('y, v, a','FontSize',14);
xlabel('t(sec)','FontSize',14);
h=legend('y','v','a',0); set(h,'FontSize',14)
