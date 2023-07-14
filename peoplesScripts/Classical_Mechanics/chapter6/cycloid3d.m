%cycloid3d.m 
%solves the motion of a charged particle in an electromagnetic 
%field in 3 dimensions. Here, x(1)=x, x(2)=vx, x(3)=y, x(4)=vy,
%x(5)=z, x(6)=vz
clear;
q=1.6e-19; m=1.67e-27; qm=q/m; %charge, mass, & ratio
t0=0.0;tmax=10.0;      %time range
%-------------------------
%Use the next 3 lines to reproduce the 2D motion case (use plot command below)
%B=[0.0;0.0;3.13e-8];  %Bx,By,Bz
%E=[5.5e-8;0.0;0.0];   %Ex,Ey,Ez
%ic=[0.0;5.0;0.0;0.0;0.0;0.0];%x0,vx0,y0,vx0,z0,vz0 - init conditions
%-------------------------
B=[1.e-8;-1.e-9;5.13e-8];    %Bx,By,Bz
E=[0.5e-8;1.e-9;-3.e-9];     %Ex,Ey,Ez
ic=[0.0;5.0;0.0;0.0;0.0;0.5];%x0,vx0,y0,vy0,z0,vz0 - init conditions
[t,r]=ode23('cycloid3d_der',[t0 tmax],ic,[],qm,B,E); %simple ode solver
%[t,r]=ode45('cycloid3d_der',[t0 tmax],ic,[],qm,B,E);% use this for a better solver
%plot(r(:,1),r(:,3),'LineWidth',1.0,'Color','black'); % 2D motion case
subplot(2,2,1)
plot3(r(:,1),r(:,3),r(:,5),'r-','LineWidth',2)%use plot3 in 3D motion
view(320,30)% can change the view angle
str=cat(2,'Charge in general E&B fields - 3D');
title(str,'FontSize',14); xlabel('x','FontSize',11); 
ylabel('y','FontSize',11); zlabel('z','FontSize',11); grid on;
axis square
subplot(2,2,2);plot(t,r(:,1));title('x vs t','FontSize',11); 
xlabel('t','FontSize',11); ylabel('x','FontSize',12); grid on;
subplot(2,2,3);plot(t,r(:,3));title('y vs t','FontSize',11); 
xlabel('t','FontSize',11); ylabel('y','FontSize',12); grid on;
subplot(2,2,4);plot(t,r(:,5));title('z vs t','FontSize',11); 
xlabel('t','FontSize',11); ylabel('z','FontSize',11); grid on;
