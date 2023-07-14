%molec.m
%program to plot the solutions of the two atom molecular potential model
%using the simple, the nonlinear approximation, and the full solutions
clear;
tmax=2; % maximum time in units of tau0
xb=3/2; %equilibrium position
xi=0.2; %initial position measured from equilibrium
ic1=[xb+xi;0.0]; % initial conditions: position, speed
%Use MATLAB's Runge-Kutta (4,5) formula
%[t,x2]=ode45('molec_der',[0.0,tmax],ic1,[]);%numeric soln (default tolerances)
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-4);   %user set Tolerances
[t,x2]=ode45('molec_der',[0.0,tmax],ic1,opt);%numerical solution
%the non0linear approximation uses A1 for amplitude
% A1 is such that at t=0 x1-xb=xi
A1=9*(-1+sqrt(1+32*xi/9))/16;
x1=xb+A1*cos(2*pi*t)+4*A1^2*(3-cos(4*pi*t))/9;
%the SHM case
x0=xb+xi*cos(2*pi*t);
plot(t,x0,'r:',t,x1,'k-.',t,x2(:,1),'b-');%position versus time
h=legend('SHO','Nonlinear Approx','Full Solution',1); set(h,'FontSize',12)
line([0,max(t)],[xb,xb],'color','black');%the equilibrium position
text(0.8,xb-0.01,'Equilibrium','FontSize',14);
str=cat(2,'x_i=',num2str(xi,3));
text(0.2,xb+xi+0.025,str,'FontSize',14);
xlabel('Time (\tau_0)','FontSize',14)
ylabel('Position (a_0)','FontSize',14)
title('Comparison of Solutions','FontSize',14)
