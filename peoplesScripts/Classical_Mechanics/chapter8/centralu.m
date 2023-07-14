%centralu.m
%program to plot the solution of a body of mass m under the action of
%a central force of the form -a*r^p where r=1/u
clear;
a=108.0;             %force strength in N/m^p
p=1;                 %power of r
r0=1; v0=6.0; hmin=0;%initial position(m), tangential velocity(m/s), angle(rad)
m=1;                 %object's mass in kg
L=m*v0*r0;           %angular momentum, v0=tangential velocity
u0=1/r0;             %initial value of u
uv0=-m*v0/L;         %initial value of the time derivative of u
hmax=2*pi;           %maximum angle
%--------------
%r(theta) as obtained from the u(theta) diff. eq.
ic1=[u0;uv0];        %initial conditions: position, initial u-dot
opt1=odeset('AbsTol',1.e-7,'RelTol',1.e-4);%user set Tolerances
[h,v]=ode45('centralu_der',[hmin,hmax],ic1,opt1,a,L,m,p);%with default tolerance
r=1./v(:,1);
%--------------
% r(theta) as obtained from the r(theta) diff. eq.
ic2=[r0;v0;hmin]; tmax=0.6;% init conds. and orbit period for tmax
[t,w]=ode45('central_der',[0.0,tmax],ic2,[],a,L,m,p);%with default tolerance
%--------------
%r(theta) comparison plots for the two methods
str=cat(2,'Central Force=-a*r^p',' with p=',num2str(p,3),', Two Methods Compared');
subplot(1,2,1); plot(h,r,'b-'); hold on
plot(w(:,3),w(:,1),'r.','MarkerSize',5)
xlabel('\theta (rad)','FontSize',14);ylabel('r(\theta)','FontSize',14)
title(str,'FontSize',14), axis([0 6 0.4 1.3])
%plot r(theta) in polar coordinates
subplot(1,2,2);polar(h,r,'b-'); hold on
m=length(w(:,3));polar(w(1:3:m,3),w(1:3:m,1),'r.');%plot every 3rd point
title('Polar Plot Comparison of Two Methods','FontSize',12)