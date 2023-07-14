%central.m
%program to plot the solution of a body of mass m under the action of
%a central force of the form -a*r^p
clear;
a=108.0;             %force strength in N/m^p
p=1;                 %power of r
tmax=1.5;            %maximum time in seconds
r0=1; v0=6.0; th0=0; %initial position(m), tangential velocity(m/s), angle(rad)
m=1;                 %object's mass in kg
L=m*v0*r0;           %angular momentum, v0=tangential velocity
ic1=[r0;v0;th0];     %initial conditions: position, v0, angle
%Use MATLAB's Runge-Kutta (4,5) formula (uncomment/comment as needed)
%opt=odeset('AbsTol',1.e-7,'RelTol',1.e-4);     %user set Tolerances
%[t,w]=ode45('central_der',[0.0,tmax],ic1,opt,a,L,m,p);%with set tolerance
[t,w]=ode45('central_der',[0.0,tmax],ic1,[],a,L,m,p);%with default tolerance
str=cat(2,'Central Force=-a*r^p',' with p=',num2str(p,3));
subplot(2,2,1)
plot(t,w(:,1),'b:',t,w(:,2),'r-')
h=legend('r','v',0); set(h,'FontSize',12)
xlabel('time (sec)','FontSize',14);
ylabel('r(t), v(t)','FontSize',14);title(str,'FontSize',14)
subplot(2,2,2)
plot(t,w(:,3))
xlabel('time (sec)','FontSize',14),ylabel('\theta (rad)','FontSize',14)
subplot(2,2,3)
plot(w(:,3),w(:,1))
xlabel('\theta (rad)','FontSize',14);ylabel('r(\theta)','FontSize',14)
subplot(2,2,4)
polar(w(:,3),w(:,1)) %plot r(theta) in polar coordinates