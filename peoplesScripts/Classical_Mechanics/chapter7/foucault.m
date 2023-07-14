%foucault.m incorporates the solution to the Foucault pendulum formulas
clear;
lth=34;             % latitude angle in degrees
th=lth*2*pi/360;    % convert to radians
g=9.8;              % acceleration due to gravity
tf=3600;            % pendulum's period is one hour for convenience
L=(tf/2/pi)^2*g;    % pendulum's length
wf=2*pi/tf;         % Foucault pendulum frequency same as sqrt(g/L)
w0=2*pi/24/3600;    % Earth's rotational frequency
w1=w0*sin(th);      % Precession frequency
w2=sqrt(w1^2+wf^2); % roughly pendulum frequency
x0=1.0; y0=0.0;     % Initial conditions
tp=abs(2*pi/w1);    % precession period
tmax=24*3600;       % 24 hour time limit
ts=tmax/1000;       % time step in seconds
t=[0:ts:tmax]; n=length(t); nh=round(n/2);
x=x0*cos(w1*t).*cos(w2*t)+y0*cos(w1*t).*sin(w2*t);
y=-x0*sin(w1*t).*cos(w2*t)+y0*cos(w1*t).*cos(w2*t);
str1=cat(2,'Latitude =',num2str(lth),' Degrees');
str2=cat(2,'\tau_p =',num2str(tp/3600),' hrs');
plot(x(1:nh),y(1:nh),'b') %first 12 hours as blue
hold on; plot(x(nh:n),y(nh:n),'r') %2nd 12 hours as red
title('24 Hour Foucault Pendulum Precession Plot','FontSize',14)
ylabel('y\prime','FontSize',14); xlabel('x\prime','FontSize',14);
text(min(x)*(1-0.05),max(y)*(1-0.1),str1,'FontSize',10,'Color','blue');
text(max(x)*(1-0.45),max(y)*(1-0.1),str2,'FontSize',10,'Color','blue');
