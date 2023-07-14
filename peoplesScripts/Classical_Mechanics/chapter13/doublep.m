%doublep.m, program to solve the double pendulum equations
%of motion numerically and plot their solutions
clear;
L1=1; L2=2; m1=1; m2=2; g=9.8;      %lengths, masses, gravity
tau=sqrt(g/(L1+L2));                %a time unit
tmax=3*tau; ts=0.05;                %simulation run time and time interval
tr=(0.0:ts:tmax);N=length(tr);      %time range, array size
th10=pi/4; th20=pi/3;               %init angle values in rad
th10d=0.0; th20d=0.0;               %init angular speeds in rad/sec
ic1=[th10;th10d;th20;th20d];        %initial conditions:
%Use MATLAB's Runge-Kutta (4,5) formula (uncomment/comment as needed)
%opt=odeset('AbsTol',1.e-8,'RelTol',1.e-5);          %user set Tolerances
%[t,w]=ode45('doublep_der',tr,ic1,opt,L1,L2,m1,m2,g);%with set tolerance
[t,w]=ode45('doublep_der',tr,ic1,[],L1,L2,m1,m2,g);  %default tolerance     
%Next: plots of the angles versus time
%w(1):theta1, w(2):theta1_dot, w(3):theta2, w(4):theta2_dot
str1=cat(2,'Double Pendulum: \theta_{10}=',num2str(th10,3),...
        ', \theta_{20}=',num2str(th20,3),' rad, d(\theta_{1}/dt)_0=',...
        num2str(th10d,3),', d(\theta_{2}/dt)_0=',num2str(th20d,3));
str2=cat(2,'L_1, L_2=',num2str(L1,3),', ',num2str(L2,3),' m');
str3=cat(2,'m_1, m_2=',num2str(m1,3),', ',num2str(m2,3),' kg',...
           ', g=',num2str(g,3),' m/s^2');
subplot(2,1,1), 
plot(t,w(:,1),'k-',t,w(:,3),'b:'), 
xlabel ('t','FontSize',13), ylabel('\theta_1, \theta_2','FontSize',13)
title(str1,'FontSize',12)
ym=min([w(:,1);w(:,3)]); yp=max([w(:,1);w(:,3)]);    %window size
axis([0,tmax,ym*(1+0.3),yp*(1+0.3)])
text(.1,yp*(1+0.05),str2)
text(.1,ym*(1+0.05),str3)
h=legend('\theta_1','\theta_2',-1); set(h,'FontSize',13)
subplot(2,1,2), plot(w(:,1),w(:,3),'r') 
xlabel ('\theta_1','FontSize',13),ylabel('\theta_2','FontSize',13)
% ================== Simulation next =================================
%Coordinates versus time
x1=L1*sin(w(:,1)); y1=L1*cos(w(:,1));
x2=x1+L2*sin(w(:,3)); y2=y1+L2*cos(w(:,3));
%r1=sqrt(x1.^2+y1.^2);r2=sqrt(x2.^2+y2.^2); %length conservation check
%figure, plot(t,[r1,r2])                    %uncomment as desired
%support is at L1+L2, where y's are measured from
v=L1+L2; y1=v-y1; y2=v-y2;  %shift the y's for simulation purpose
figure, vx=max([v;x1;x2]);vy=max([v;y1;y2]);
axis([-vx,vx,0,vy])                                  %window size
 for i=1:N
     clf
     axis([-vx,vx,0,vy])
     hold on
     plot(0,v,'ko');                                 %pivot point
     line([0,x1(i)],[v,y1(i)],'color', 'k');         %arm1
     h(1)=plot(x1(i),y1(i),'k.');                    %m1 at x1,y1
     line([x1(i),x2(i)],[y1(i),y2(i)],'color','b');  %arm2 
     h(2)=plot(x2(i),y2(i),'b.');                    %m2 at x2,y2
     pause(.05)
 end
h(3)=plot(x1,y1,'k:');                               %m1 trace
h(4)=plot(x2,y2,'b-.');                              %m2 trace
h=legend(h,'m_1','m_2','m_1 trace', 'm_2 trace');
set(h,'FontSize',13)
xlabel('x_1, x_2','FontSize',13), ylabel('y_1, y_2','FontSize',13)
title(str1,'FontSize',12), 
text(-v*(1-0.1),v*(1-.05),str2)
text(-v*(1-0.1),v*(1-.1),str3)