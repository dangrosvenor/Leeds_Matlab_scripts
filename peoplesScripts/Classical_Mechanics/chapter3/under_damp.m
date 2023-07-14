%under_damp.m
%plots the underdamped HO solution
clear;
m=0.05;                             %mass
k=1;                                %spring constant
c=0.08;                             %drag coefficient
x0=1.0;                             %initial position
v0=5.0;                             %initial speed
gam=c/2/m;                          %critical gamma
wo=sqrt(k/m);                       %SHO natural frequency
desc=wo^2-gam^2;                    %must be positive this time
if desc <= 0;                       %ensure appropriate problem conditions
    disp('gam needs to be smaller'); 
    return; 
end
w=sqrt(desc);                       %underdamped frequency
B=sqrt(x0^2+(v0+gam*x0)^2/w^2);     %constant B for underdamped
th=atan(w*x0/(v0+gam*x0));          %angle theta
tmax=5;                             %maximum time
NPTS=100;                           %number of points
t=[0:tmax/NPTS:tmax];               %time array
x=B*exp(-gam*t).*sin(w*t+th);       %underdamped solution
xe=B*exp(-gam*t);                   %the decay envelope
plot(t,x,'b:',t,xe,'r-.');
title('Underdamped Harmonic Oscillator','FontSize',14)
ylabel('Underdamped, Decay Envelope','FontSize',14);
xlabel('t','FontSize',14);
h=legend('Underdamped','Decay Envelope',1); set(h,'FontSize',14);