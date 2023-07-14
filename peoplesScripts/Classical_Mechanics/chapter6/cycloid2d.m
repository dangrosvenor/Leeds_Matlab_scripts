%cycloid2d.m - plots the motion of a charged particle in the presence of
%electric (+x direction) and magnetic (+z direction) fields.
clear;
tmax=10;
q=1.6e-19; m=1.67e-27; %proton charge(Coulombs), proton mass (kg)
B0=3.13e-8; E0=5.5e-8; %B and E fields in Tesla & volts/meter
v0x=5; v0y=0.0; %initial x,y speeds in m/s
x0=0.0; y0=0.0; %initial x,y position (m)
f=atan(v0x/(v0y+E0/B0)); %the angle fi
wc=q*B0/m;               %cyclotron frequency
A=-sqrt(v0x^2+(v0y+E0/B0)^2)/wc;%(pick - sign for convenience) constant A
a=x0-A*cos(f); b=y0+A*sin(f);   %constants based on init conditions
t=[0:0.01:tmax];
x=A*cos(wc*t+f)+a; 
y=-A*sin(wc*t+f)-E0*t/B0+b;
plot(x,y,'b-.','LineWidth',1.5)
str=cat(2,'Motion of a charge in E(+x-east), and B(+z-out of paper) fields');
title(str,'FontSize',13), ylabel('y','FontSize',14)
xlabel('x','FontSize',14); grid on;
