%pend2.m
%program to plot the solutions of the harmonic oscillator using the simple, the
%nonlinear approximation, and the full sin(theta) term solution
clear;
cf=2*pi/360; %conversion factor from degrees to radians
w0=1; % let the frequency of the SHM be one
tau0=2*pi/w0; %period for the SHM
tmax=4*tau0; % maximum time
th=90; %initial amplitude in degrees
thr=th*cf;% initial angle in radians
ic1=[thr;0.0]; % initial conditions for angle and angular speed
%Use MATLAB's Runge-Kutta (4,5) formula
%[t,th2]=ode45('pend2_der',[0.0,tmax],ic1,[],w0);%numerical solution (default tolerances)
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-4);      %user set Tolerances
[t,th2]=ode45('pend2_der',[0.0,tmax],ic1,opt,w0);%numerical solution
%the nonlinear approximation
om=w0*sqrt(1-thr^2/8);
a3=w0^2/6;
A1=thr;
A3=a3*A1^3/(27*a3*A1^2-32*w0^2);
th1=thr*cos(om*t)+A3*cos(3*om*t);
%the SHO case
th0=thr*cos(w0*t);
plot(t,th0/cf,'r:',t,th1/cf,'k-.',t,th2(:,1)/cf,'b-');%Amplitude in degrees
h=legend('Standard SHO','Nonlinear Approx','Full Solution',1);
set(h,'FontSize',12)
axis([0 max(t) -th th*(1+0.4)]);
str=cat(2,'\theta_0=',num2str(th,3));
text(5,th*(1+0.2),str,'FontSize',12);
xlabel('Time (sec)','FontSize',14)
ylabel('Amplitude (degrees)','FontSize',14)
title('Comparison of Solutions','FontSize',14)
