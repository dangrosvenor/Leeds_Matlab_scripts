%pend1.m
%Program designed to compare the period of a pendulum with the next approximation
%beyond the linear simple pendulum formula. In the simple pendulum, the period is 
%independent of amplitude
clear;
N=100;
thmax=90;% note at 180, the period is infinite - pendulum gets stuck
dth=thmax/N;
th=[0:dth:thmax];% in degrees
%MATLAB defines the complete elliptic integral of the first kind as:
% ellipke(m)=Integral of (1-m*sin^2(t)) on the interval 0 <= t <= pi/2
% Thus Matlab's m is actually our k^2, where k=sin(theta/2)
m=sin(th*2*pi/360/2).^2 ;
% Periods in units of simple pendulum tau=2*pi*sqrt(l/g)
y1=1./sqrt(1-(th*2*pi/360).^2/8);    %the nonlinear approximation
y2=2*ellipke(m)/pi;                  %full solution for the period versus amplitude
plot(th,y1,'b.',th,y2,'k-')
line([0,thmax],[1,1],'color','red'); %the simple harmonic oscilator period taken as one
axis([0 thmax 0.95 1.2]);
h=legend('\tau non-linear-approx','\tau-full','\tau_0',0);
set(h,'FontSize',14)
xlabel('Amplitude(Degrees)','FontSize',14)
ylabel('Period','FontSize',14)
title('Comparison of the Periods','FontSize',14)
