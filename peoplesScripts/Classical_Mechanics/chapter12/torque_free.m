%torque_free.m - plots the angular frequency and momentum for torque free
%motion of a top versus time in the body (S') frame
clear; ap=180/pi;                             %angle conversion
I1=1.3; I2=I1; I3=1.5; gam=(I3-I1)/I1;        %principal moments of inertia
w3=1;fi=0; c=1;                               %initial values, c, fi, w3
omb=gam*w3;                                   %precessional freq.
tmax=2*pi/omb; N=40; ts=tmax/N;               %time limits
t=[0:ts:tmax];                                %time range
w1=c*cos(omb*t+fi);w2=c*sin(omb*t+fi);        %w1, w2 versus time                  
L1=I1*w1; L2=I2*w2; L3=I3*w3;                 %Ang. Mom. components
v1=max(L1); v2=max(L2); v3=max(L3*(1+0.1));   %view window
fi_b=atan(c/w3);                              %angle between w and w3
fi_L=acos((I1*c^2+I3*w3^2)/...                %angle between w and L
         (sqrt((c^2+w3^2)*(I1^2*c^2+I3^2*w3^2))));
fi_s=atan(I1*tan(fi_b)/I3);                   %angle between L and w3
omL=omb*sin(fi_b)/sin(fi_L);                  %prec. freq. of w about L
%omL=sqrt(I1^2*c^2+I3^2*w3^2)/I1              %another formula for omL
w=sqrt(c^2+w3^2);                             %magnitude of w
L=sqrt(I1^2*c^2+I3^2*w3^2);                   %magnitude of L
for i=1:N
    clf
    axis ([-v1,v1,-v2,v2,0,v3])
    h(1)=line([0,w1(i)],[0,w2(i)],[0,w3],'color', 'k',...
             'LineStyle','-.','linewidth', 1.5);                %w  line
    h(2)=line([0,L1(i)],[0,L2(i)],[0,L3],'color', 'r',...
             'LineStyle','--','linewidth', 1.5);                %L  line
    h(3)=line([0,0],[0,0],[0,w3],'color', 'b', 'linewidth',1.5);%w3 line
    box on
    pause(.1)
end
hold on
h(4)=plot3 (w1,w2,w3*(t+0.01)./(t+0.01),'k-.'); %plot the w cone (body)
h(5)=plot3 (L1,L2,L3*(t+0.01)./(t+0.01),'r:');  %plot the L cone (space)
hh=legend(h,'\omega','L','\omega_3','Body Cone','Space Cone',4);
set(hh,'FontSize',8,'Position',[0.7 0.3 0.2 0.17])
str1=cat(2,'Torque Free Motion of a Top, I_1=',num2str(I1,3),...
           'kgm^2, I_3=',num2str(I3,3),'kgm^2',', \omega=',num2str(w,3)...
           ,'rad/s, L=',num2str(L,3),'kgm^2/s');
title(str1,'FontSize',10), xlabel('x','FontSize',14)
ylabel('y','FontSize',14), zlabel('z','FontSize',14)
str2=cat(2,'\phi_b=',num2str(fi_b*ap,3),'^o, \phi_L= ',...
            num2str(fi_L*ap,3),'^o, \phi_s=',num2str(fi_s*ap,3),'^o');
str3=cat(2,'\Omega_b=',num2str(omb,3),'rad/s, \Omega_L=',...
            num2str(omL,3),'rad/s');
text (-v1*(1-0.1),v2,v3*(1-0.75),str2,'FontSize',11)
text (-v1*(1-0.1),v2,v3*(1-0.9),str3,'FontSize',11)
