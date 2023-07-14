%fixed_axis.m -animates the position of a rod-mass system, its angular momentum,
%and the related torque
clear; warning off;
w=pi;tau=2*pi/w;      %angular velocity and period
t=[0:tau/20:1.5*tau]; %time variable
N=length(t);          %number of points
th=25;                %angle of rod tilt from z axis
a=1; thr=th*pi/180;   %rod length, angle in radians
s=sin(thr); c=cos(thr);R=a*s; z=a*c; co=cot(thr);%for orbit radius and height
m=1;I=m*R^2;A=I*w;    %mass, moment if inertia, angular momentum amplitude
v=max([R,z,A]);       %view window parameter
vxy=max([R,z,A*w*co]);%view window parameter
for i=1:N
clf
axis ([-vxy,vxy,-vxy,vxy,0,v])
ct=cos(w*t(i)); st=sin(w*t(i)); %needed to get x,y coords vs time
x=R*ct; y=R*st;                 %x,y particle coordinates versus time
line([0,x],[0,y],[0,z],'color', 'black', 'linewidth', 1)%massless rod line
line([x,x],[y,y],[z,z],'color', 'black','LineStyle','.','linewidth', 1,...
    'Marker','.','MarkerSize',20)       %the particle
Lx=-A*co*ct; Ly=-A*co*st; Lz=A;         %angular momentum components
Tx=A*w*co*st; Ty=-A*w*co*ct; Tz=0.0;    %torque components
line([0,Lx],[0,Ly],[0,Lz],'color','red','linewidth',1.5,...
    'LineStyle','--') %L (ang. Mom.)
line([0,Tx],[0,Ty],[0,Tz],'color','green', 'linewidth',1.5,...
    'LineStyle','-.')%torque
line([0,0],[0,0],[0,z],'color','m','LineStyle',':','linewidth', 1.5)%z axis
line([Lx,Lx],[Ly,Ly],[Lz,Lz],'color','red','Marker','d','MarkerSize',5)  %Arrow for L
line([Tx,Tx],[Ty,Ty],[Tz,Tz],'color','green','Marker','d','MarkerSize',5)%Arrow for T
grid on
pause(.1)
end
h=legend('rod','mass','L','\tau','z-axis',2); set(h,'FontSize',14)
hold on
plot3(R*cos(w*t),R*sin(w*t),z*t./t,':')% particle trace
str=cat(2,'Fixed-axis massless rod-point mass rotating system',...
          ', \theta = ',num2str(th,3));
title(str,'FontSize',14)
xlabel('x','FontSize',14), ylabel('y','FontSize',14), zlabel('z','FontSize',14)
%r=[x,y,z], L=[Lx,Ly,Lz], dot(r,L) %use to check if r, L are perpendicular