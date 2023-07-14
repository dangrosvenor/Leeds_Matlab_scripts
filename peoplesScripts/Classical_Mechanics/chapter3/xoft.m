%xoft.m
% plots x=r*cos(2pi*t/tau)
clear;
tau=2;               %period
NPTS=100;            %number of points
t=[0:2/NPTS:2*tau];  %time array
r=2.0;               %amplitude
x=r*cos(2*pi*t/tau); %position array
plot(t,x);
title('Plot of x(t)=r\cdotcos(2\pit/\tau) vs time','FontSize',14)
ylabel('x(t)','FontSize',14);
xlabel('t(sec)','FontSize',14);
