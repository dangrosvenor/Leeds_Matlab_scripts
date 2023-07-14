%molec_mu.m
%program to plot the coordinates of the atoms of a molecule. The molecular 
%potential model used is a simple 1/x^4-1/x^3 type. The method of successive 
%approximations has been used for the relative coordinate. The center of 
%mass is under the action due to gravity.
clear;
xcm=0.0;              %x center of mass position
v0y=0.0;y0=20.0;g=9.8;%initial y-velocity, and y-position, g=gravity
tmax=v0y/g+sqrt((v0y/g)^2+2*y0/g);%time to reach ground is used for tmax
xb=3/2;         %equilibrium relative coordinate (molecule bond length)
xi=0.2;         %initial position measured from equilibrium
t=[0:0.01:tmax];%time array
A1=9*(-1+sqrt(1+32*xi/9))/16;                 %A1, at t=0 x-xb=xi
x=xb+A1*cos(2*pi*t)+4*A1^2*(3-cos(4*pi*t))/9; %relative coordinate
ycm=y0+v0y*t-0.5*g*t.^2;                      %y-center of mass position
x1=xcm-0.5*x; x2=xcm+0.5*x;                   %coordinates of the atoms
plot (x1,ycm,'r-.')                           %1st atom x,y positions
hold on
plot (x2,ycm,'k')                             %2nd atom x,y positions
xlabel('x (a_0)','FontSize',13), ylabel('y (a_0)','FontSize',13)
title('Plot of atomic x positions as the molecule free falls','FontSize',14)
line([-xb/2,-xb/2],[0,y0],'color','blue','LineStyle',':') %atom1 ave x
line([xb/2,xb/2],[0,y0],'color','blue','LineStyle',':')   %atom2 ave x
line([-xb/2,-xb/2*(1-0.6)],[y0/2,y0/2],'LineStyle',':')   %bond length
line([xb/2,xb/2*(1-0.6)],[y0/2,y0/2],'LineStyle',':')     %bond length
str=cat(2,' bond length=',num2str(xb,3),'a_0');
text(-xb/2*(1-0.6),y0/2,str,'FontSize',14);
h=legend('1st atom','2nd atom',0); set(h,'FontSize',14)
