%torque_free_s.m - plots the angular frequency and momentum for 
%torque free motion of a top versus time in the body (S') as well 
%as in the space frame (S)
clear; ap=180/pi;                              %angle conversion
I1=1.3; I2=I1; I3=1.8; gam=(I3-I1)/I1;         %principal moments of inertia
w3=1; fi=0; c=1;                               %initial values, c, fi, w3
fi_b=atan(c/w3);                               %angle between w and w3
fi_L=acos((I1*c^2+I3*w3^2)/...                 %angle between w and L
         (sqrt((c^2+w3^2)*(I1^2*c^2+I3^2*w3^2))));
fi_s=atan(I1*tan(fi_b)/I3);                    %angle between L and w3
omb=gam*w3;                                    %precessional freq.
omL=omb*sin(fi_b)/sin(fi_L);                   %prec. freq. of w about L 
tmax=.9*2*pi/omb; N=100; ts=tmax/N;            %time limits
t=[0:ts:tmax];                                 %time range
w1=c*cos(omb*t+fi);w2=c*sin(omb*t+fi);         %w1, w2 versus time
L1=I1*w1; L2=I2*w2; L3=I3*w3;                  %Ang. Mom. components
w=sqrt(c^2+w3^2);                              %magnitude of w
L=sqrt(I1^2*c^2+I3^2*w3^2);                    %magnitude of L
%=== shift old e3 by fi_s to swap places with old L, preserve magnitude
e31=1*sin(fi_s)*cos(omL*t+fi);                 %must move with omL
e32=1*sin(fi_s)*sin(omL*t+fi);                 %must move with omL
e33=1*cos(fi_s);
%new w and components, shift by fi_s, preserve magnitude, change rate
c=w*sin(fi_b-fi_s);
w3=w*cos(fi_b-fi_s);
w1=c*cos(omL*t+fi);w2=c*sin(omL*t+fi);         %use omL precession rate
%new L now in old w3 direction and fixed, preserve magnitude
L1=0; L2=0; L3=L;
v1=max(e31); v2=max(e32); v3=max(L3*(1+0.1));  %view window
for i=1:N
    clf
    axis([-v1,v1,-v2,v2,0,v3])
    hold on
    h(1)=line([0,L1],[0,L2],[0,L3],'color', 'r',...        %L  line
             'LineStyle','--','linewidth', 1.5); 
    h(2)=line([0,w1(i)],[0,w2(i)],[0,w3],'color','m',...   %S frame w line     
        'LineStyle','--','linewidth', 1);
    h(3)=line([0,e31(i)],[0,e32(i)],[0,e33],'color','b',...%e3 line
              'linewidth',1.5);
    box on
    pause(.05)
end
h(4)=plot3 (w1,w2,w3*(t+0.01)./(t+0.01),'m:');     %plot the S w cone
hh=legend(h,'L','\omega S-frame','e_3','\omega cone',1);
str1=cat(2,'Torque Free Motion of a Top, I_1=',num2str(I1,3),...
           'kgm^2, I_3=',num2str(I3,3),'kgm^2');
str1=cat(2,'Torque Free Motion of a Top, I_1=',num2str(I1,3),...
         'kgm^2, I_3=',num2str(I3,3),'kgm^2',', \omega=',num2str(w,3),...
         'rad/s, L=',num2str(L,3),'kgm^2/s');       
title(str1,'FontSize',10), xlabel('x','FontSize',14)
ylabel('y','FontSize',14), zlabel('z','FontSize',14)
str2=cat(2,'\phi_b=',num2str(fi_b*ap,3),'^o, \phi_L= ',...
            num2str(fi_L*ap,3),'^o, \phi_s=',num2str(fi_s*ap,3),'^o');
str3=cat(2,'\Omega_b=',num2str(omb,3),'rad/s, \Omega_L=',...
            num2str(omL,3),'rad/s');
text (-v1*(1-0.4),v2,v3*(1-0.15),str2,'FontSize',11)
text (-v1*(1-0.4),v2,v3*(1-0.05),str3,'FontSize',11)
