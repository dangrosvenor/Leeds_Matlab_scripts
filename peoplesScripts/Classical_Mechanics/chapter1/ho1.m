%ho1.m
%Calculation of position, velocity, and acceleration for a harmonic
%oscillator versus time. The equations of motion are used for small time intervals
clear;
%NPTS=100;TMAX=1.0;%example Maximum number of points and maximum time
TTL=input(' Enter the title name TTL:','s');%string input
NPTS=input(' Enter the number calculation steps desired NPTS: ');
TMAX=input(' Enter the run time TMAX: ');
NT=NPTS/10;%to print only every NT steps
%K=1000;M=5.0;C=0.0;E=0.0;W=0.0;x0=0.1;v0=0.0;% example Parameters
K=input(' Enter the Spring contant K: ');
M=input(' Enter the bob mass M: ');
C=input(' Enter the damping coefficient C: ');
E=input(' Enter the magnitude of the driving force E: ');
W=input(' Enter the driving force frequency W: ');
x0=input(' Enter the initial position x0: ');% Initial Conditions
v0=input(' Enter the initial velocity v0: ');% Initial Conditions
t0=0.0;% start at time t=0
dt=TMAX/NPTS;%time step size
fprintf(' Time step used dt=TMAX/NPTS=%7.4f\n',dt);%the time step being used
F=-K*x0-C*v0+E*sin(W*t0); % initial force
a0=F/M;% initial acceleration
fprintf('    t       x       v       a\n');%output column labels
v(1)=v0;
x(1)=x0;
a(1)=a0;
t(1)=t0;
fprintf('%7.4f %7.4f %7.4f %7.4f\n',t(1),x(1),v(1),a(1));%print initial values
for i=1:NPTS
    v(i+1)=v(i)+a(i)*dt;                   %new velocity
    x(i+1)=x(i)+v(i+1)*dt;                 %new position
    t(i+1)=t(i)+dt;                        %new time
    F=-K*x(i+1)-C*v(i+1)+E*sin(W*t(i+1));  %new force
    a(i+1)=F/M;                            %new acceleration
% print only every NT steps
    if(mod(i,NT)==0)
        fprintf('%7.4f %7.4f %7.4f %7.4f\n',t(i+1),x(i+1),v(i+1),a(i+1));
    end;
end;
subplot(3,1,1)
plot(t,x,'k-');
ylabel('x(t) (m)','FontSize',14);
h=legend('position vs time'); set(h,'FontSize',14);
title(TTL,'FontSize',14);
subplot(3,1,2)
plot(t,v,'b-');
ylabel('v(t) (m/s)','FontSize',14);
h=legend('velocity vs time'); set(h,'FontSize',14)
subplot(3,1,3)
plot(t,a,'r-');
ylabel('a(t) (m/s^2)','FontSize',14);
xlabel('time (sec)','FontSize',14);
h=legend('acceleration vs time'); set(h,'FontSize',14)