%over_crit_damp.m
%plots the overdamped and critically damped HO solutions
clear;
m=0.05;                             %mass
k=1;                                %spring constant
c=0.5;                              %drag coefficient
x0=1.0;                             %initial position
v0=5.0;                             %initial speed
gam=c/2/m;                          %critical gamma
desc=gam^2-k/m;                     %must be positive
if desc <= 0;                       %ensure appropriate problem conditions
    disp('gam needs to be greater'); 
    return; 
end
gam1=gam+sqrt(desc);                %overdapmed gamma1
gam2=gam-sqrt(desc);                %overdapmed gamma2
Bo=(v0+gam1*x0)/(gam1-gam2);        %constant B for overdamped
Ao=x0-Bo;                           %constant A for overdamped
Ac=x0;                              %constant A for critically damped 
Bc=v0+gam*x0;                       %constant B for critically damped
tmax=2;                             %maximum time
NPTS=100;                           %number of points
t=[0:tmax/NPTS:tmax];               %time array
xo=Ao*exp(-gam1*t)+Bo*exp(-gam2*t); %overdamped
xc=Ac*exp(-gam*t)+Bc*t.*exp(-gam*t);%critically damped
plot(t,xo,'b:',t,xc,'r-.');
title('Overdamped and Critically Damped Comparison','FontSize',14)
ylabel('Xoverdamped, Xcritically-damped','FontSize',14);
xlabel('t','FontSize',14);
h=legend('Overdamped','Critically Damped',1); set(h,'FontSize',14);