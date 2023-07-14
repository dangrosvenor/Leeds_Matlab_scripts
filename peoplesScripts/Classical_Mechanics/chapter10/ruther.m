%ruther.m - program to simulate the Rutherford Scattering alpha particle trajectory
%numerically, given the initial projectile energy and impact parameter.
clear;
rcm=0.0;                 %target position
za=2; zt=79;             %za=projectile, zt=target charges (2->alpha, 79->gold)
K=za*zt;                 %dimensionless force constant
m=1;                     %projectile mass in units of alpha particle mass
vb=0.01965;              %velocity units (in inits of c=light speed)
ma=3730e6;               %alpha particle mass energy in eV
Ene=5e6;                 %initial projectile energy in eV
v0=sqrt(2*Ene/m/ma)/vb;  %initial speed in units of vb
b=20;                    %impact parameter
L=100;                   %for window viewing purpose
x10=60*L;                %initial x position - far away (a_b=1 Fermi units)
y10=b;                   %initial y position is impact parameter
v1x0=-sign(x10)*sqrt(v0^2-K/x10/m/2); %initial x speed based on energy conservation
                                      %toward target, opposite sign to x10
v1y0=0.0;                %zero y-speed - projectile coming in horizontally
tmax=2500;               %max calculation time in units of tau_b (tau_b~1.7e-22 sec)
ic1=[x10,v1x0,y10,v1y0]; %initial conditions
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-6);%user set Tolerances
[t,w]=ode45('ruther_der',[0.0,tmax],ic1,opt,K,m);%with set tolerance
% str=cat(2,'y versus x - binary system tau=',num2str(tau,4),' tau_0');
[nr,nc]=size(w);%get rows and columns of r
%for i=1:nr
 for i=round(7*nr/16):round(13*nr/16) %not all points animated
    clf;
    hold on
    plot(rcm,rcm,'*')
    plot(w(:,1),w(:,3),'k:') %plot y versus x
    axis ([-L/2 L -L/2 L])
    plot(w(i,1),w(i,3),'ro','MarkerSize',10)
    pause(0.0125)
 end
xlabel('x (a_b)','FontSize',13), ylabel('y (a_b)','FontSize',13)
str=cat(2,'E_k= ',num2str(Ene,3),' eV',', b= ',num2str(b,3),' a_b',...
    ', m= ',num2str(m,3),' m_\alpha');
text(-L/2*(1-0.1),L*(1-0.1),str,'FontSize',10)
%== asymptotes next ==
x=w(:,1);y=w(:,3);                  %let x, y take their names
lx=length(x);                       %x array length
sa=(y(lx)-y(lx-1))/(x(lx)-x(lx-1)); %outgoing asymptote slope (use last point)
xa=-L/2:L;                          %variable to plot outgoing asymptote
ya=y(lx)+sa*(xa-x(lx));             %asymptote line
plot (xa,ya,'b-.')                  %outgoing asymptote
line([-L/2,L],[b,b],'Color','b','LineStyle','-.') %incoming asymptote
if(sa>=0) THETA=pi-atan(sa); end    %first quadrant case
if(sa <0) THETA=pi-(atan(abs(1/sa))+pi/2); end %2nd quadrant case
r=sqrt(x.^2+y.^2); rmin=min(r);     %get r(t) and rmin
%====================
ye=-sign(K)*sqrt(1+m^2*v1x0^4*b^2/K^2);% eccentricity
str2=cat(2,'Numerical Simulation, r_{min}= ',num2str(rmin,3),' a_b, \Theta=',...
    num2str(THETA,3),' rad = ',num2str(THETA*360/2/pi,3),'^o',', e= ',num2str(ye,3));
title(str2,'FontSize',11)
line([-L/2,L],[0,0],'Color','red','LineStyle','--')
h=legend('Gold Target','Projectile Path','Projectile Particle',' Asymptotes',4);
set(h,'FontSize',13)
