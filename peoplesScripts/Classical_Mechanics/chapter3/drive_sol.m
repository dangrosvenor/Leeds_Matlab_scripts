%drive_sol.m
%plots the solution for a driven HO and the applied external force
clear;
w=input(' Enter the driving frequency: ');
c=input(' Enter the drag coefficient value: ');
x0=1.0;                            %initial position for homogeneous part
v0=5.0;                            %initial speed for homogeneous part
m=0.5;                             %mass
k=0.5;                             %spring constant
F0=0.5;                            %driving force amplitude
theta=0;                           %driving force initial phase angle
wo=sqrt(k/m);                      %SHO natural frequency
dt=0.05;                           %time step
tau=2*pi/w;                        %force's period of rotation
tmax=5*tau;                        %maximum time in terms of tau
NPTS=tmax/dt;                      %number of points
gam=c/2/m;                         %find gamma
desc=(2*gam*w).^2+(wo^2-w^2).^2; 
A=F0/m/sqrt(desc);                 %The driven ho amplitude
den=wo^2-w^2;
if den==0, den=1.e-3; end
if(w <= wo)
  ph=atan(2*gam*w/den);         %phase difference between force and soln
else
  ph=pi+atan(2*gam*w/den);      %shift by pi needed if w > wo
end
%fprintf('gamma =%7.4f\n',gam); %uncomment to print value to screen
delta=theta-ph;                 %the forced solution's phase 
t=[0:tmax/NPTS:tmax];
xf=A*cos(w*t+delta);            %the forced or particular solution
F=F0*cos(w*t+theta);
%========================================================================
% forced solution
subplot(2,1,1)       %setup 2 x 1 matrix plot - 1st window
plot(t,F,'k-',t,xf,'b:')
                     %num2str(c,p) converts c to a string with p digits
                     %cat(2,'a','b') concatenates a and b
str=cat(2,'\gamma=',num2str(gam,3),', \omega_D=',num2str(w,3),...
    ', \phi=',num2str(ph,3));
top=max(A,F0)*(1+0.3); %to create a windows large enough
topt=max(A,F0)*(1+0.1);%to post parameters on the plot
axis([0 tmax -top  top]);
text(0,topt,str,'FontSize',12,'Color','red');
title('Driving Force and Driven HO Solution','FontSize',8)
ylabel('x_f(t), F(t)');
xlabel('t','FontSize',8);
legend('F(t)','x(t)');
%========================================================================
% homogeneous + forced solution
subplot(2,1,2)                      %setup 2 x 1 matrix plot - 2nd window
desc=wo^2-gam^2;                    %must be positive
if desc <= 0;                       %ensure homogeneous problem conditions
    disp('gam needs to be smaller'); 
    return; 
end
wu=sqrt(desc);                       %underdamped homogeneous soln frequency
th=atan(wu*x0/(v0+gam*x0));          %angle theta
B=sqrt(x0^2+(v0+gam*x0)^2/wu^2);     %constant B for underdamped
xh=B*exp(-gam*t).*sin(wu*t+th);      %homogeneous solution
plot(t,xh+xf,'r');                   %plot homogeneous+particular solution
title('Homogenous+Particular solutions','FontSize',8)
ylabel('x(t)=x_h(t)+x_f(t)');
xlabel('t','FontSize',8);
axis([0 max(t) min(xh+xf) max(xh+xf)])