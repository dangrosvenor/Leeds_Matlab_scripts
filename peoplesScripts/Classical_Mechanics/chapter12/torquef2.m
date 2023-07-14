%torquef2.m
%program to solve Euler's equations for an ellipsoid without
%torques
clear;
I1=145/46; I2=1589/251; I3=1242/151;                %calculated by ellipso.m
gam1=(I3-I2)/I1; gam2=(I1-I3)/I2; gam3=(I2-I1)/I3;  %gammas
tmax=12;ts=0.05;             %simulation run time and tome interval
tr=[0.0:ts:tmax];            %time range
w10=1; w20=1; w30=1;         %initial values
ic1=[w10;w20;w30;];          %initial conditions:
%Use MATLAB's Runge-Kutta (4,5) formula (uncomment/comment as needed)
%opt=odeset('AbsTol',1.e-8,'RelTol',1.e-5);           %user set Tolerances
%[t,w]=ode45('torquef2_der',tr,ic1,opt,gam1,gam2,gam3);%with set tolerance
[t,w]=ode45('torquef2_der',tr,ic1,[],gam1,gam2,gam3);  %default tolerance     
subplot(2,2,1), plot(t,w(:,1),'k')
xlabel('t','FontSize',14),ylabel('\omega_1','FontSize',14)
title(['Ellipsoid System: I_1=',num2str(I2,2),', I_2='...
        ,num2str(I2,2),', I_3=',num2str(I3,2)],'FontSize',12)
subplot(2,2,2), plot(t,w(:,2),'b')
xlabel('t','FontSize',14),ylabel('\omega_2','FontSize',14)
Title(['\gamma_1=',num2str(gam1,3),', \gamma_2='...
        ,num2str(gam2,3),', \gamma_3=',num2str(gam3,3)],'FontSize',12)
subplot(2,2,3), plot(t,w(:,3),'r')
xlabel('t','FontSize',14),ylabel('\omega_3','FontSize',14)
subplot(2,2,4), plot3(w(:,1),w(:,2),w(:,3),'color',[0.8 0.4 0.1])
view([1 1/2 1/2])          %viewpoint in x,y,z cartesian coords
v1x=min(w(:,1));v1y=min(w(:,2));v1z=min(w(:,3)); %for window view
v2x=max(w(:,1));v2y=max(w(:,2));v2z=max(w(:,3));
axis ([v1x,v2x,v1y,v2y,v1z,v2z])                 %window size
grid on; 
%box on
xlabel('\omega_1','FontSize',14),ylabel('\omega_2','FontSize',14)
zlabel('\omega_3','FontSize',14)
%================== Simulation next =================================
figure
Iv=(I1+I2+I3)/3;                                    %Average Moment
L1=I1*w(:,1)/Iv; L2=I2*w(:,2)/Iv; L3=I3*w(:,3)/Iv;  %Ang. Mom./Iv
v1a=max(L1); v2a=max(L2); v3a=max(L3);
v1b=max(w(:,1)); v2b=max(w(:,2)); v3b=max(w(:,3));
v1=max(v1a,v1b);v2=max(v2a,v2b);v3=max(v3a,v3b);    %view window
N=length(tr);
for i=1:5:N
    clf
    axis([-v1,v1,-v2,v2,0,v3])
    hold on
    h(1)=line([0,L1(i)],[0,L2(i)],[0,L3(i)],'color', 'k',...  %L/Iv line
             'LineStyle','-','linewidth', 2); 
    h(2)=line([0,w(i,1)],[0,w(i,2)],[0,w(i,3)],'color','b',...%w line     
        'LineStyle','-.','linewidth', 2);
    box on
    pause(.05)
end
h(3)=plot3(L1,L2,L3,'k:');                 %plot the L/Iv trace
h(4)=plot3(w(:,1),w(:,2),w(:,3),'b:');     %plot the w trace
hh=legend(h,'L/I_{avg}','\omega','L trace','\omega trace',1);
xlabel('x','FontSize',14),ylabel('y','FontSize',14)
zlabel('z','FontSize',14)
title(['Ellipsoid System: \gamma_1=',num2str(gam1,3),', \gamma_2='...
        ,num2str(gam2,3),', \gamma_3=',num2str(gam3,3)],'FontSize',13)
