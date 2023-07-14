%fofx.m
clear;
m=1.0;                          %mass
k=0.01;                         %spring constant
w=sqrt(k/m);                    %natural frequency
x0=0.0;                         %initial position
v0=0.5;                         %initial velocity
t=[0:0.05:2*pi/w];              %time array from zero to one oscillation period
E0=0.5*m*v0^2;                  %total initial energy
x=sqrt(2*E0/k)*sin(w.*t);       %position versus time array
v=v0*cos(w.*t);                 %velocity versus time arry
%a=-k*x/m;                      %acceleration versus time array if needed
PE=0.5*k*x.^2;                  %potential energy array
KE=0.5*m*v.^2;                  %kinetic energy array
E=PE+KE;                        %total energy array
F=-k*x;                         %force array
plot(x,PE,'k-',x,KE,'b:',x,E,'r-.',x,F,'m--');
title('Spring-Mass Simple Harmonic Energy-Force Relation','FontSize',14)
ylabel('PE, KE, E, F','FontSize',14);
xlabel('x(m)','FontSize',14);
h=legend('PE','KE','E','F',3); set(h,'FontSize',14)