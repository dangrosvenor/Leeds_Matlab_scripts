%inter_spr2.m
%plots the coordinate solutions for the full bimodal coupled spring-mass system
clear;
m=1;                            %masses
k0=1.0; k=10.0;                 %spring constants
x10=1;  x20=0;                  %initial positions
xs=(x10+x20)/2;
xd=(x10-x20)/2;
om1=sqrt(k/m);                  %frequencies
om2=sqrt((k+2*k0)/m);
om=min(om1,om2);                %lowest frequency: for time range purposes
tau=2*pi/om; tmax=20*tau;
str=cat(2,'m=',num2str(m),', k_0=',num2str(k0),', k=',num2str(k),...
', x_{10}=',num2str(x10),', x_{20}=',num2str(x20),...
', \omega_1=',num2str(om1,2),', \omega_2=',num2str(om2,2));
t=[0:tau/50:tmax];              %plotting time range
x1=xs*cos(om1*t)+xd*cos(om2*t); %solutions
x2=xs*cos(om1*t)-xd*cos(om2*t);
% -- can also write x1 and x2 as x3 and x4 if desired
%x3=x10*cos((om1+om2)*t/2).*cos((om2-om1)*t/2)...
%    +x20*sin((om1+om2)*t/2).*sin((om2-om1)*t/2);
%x4=x10*sin((om1+om2)*t/2).*sin((om2-om1)*t/2)...
%    +x20*cos((om1+om2)*t/2).*cos((om2-om1)*t/2);
%plot(t,x1,t,x2);
plot(t,x1,'b-',t,x2,'r:');
axis([0 tmax -1.2  1.2]);
text(0.5,1.05,str,'FontSize',12,'Color','black');
title('Coupled Spring-Mass Bimodal System','FontSize',14)
h=legend('x1','x2',4); set(h,'FontSize',14)
ylabel('Position','FontSize',14);
xlabel('Time','FontSize',14);