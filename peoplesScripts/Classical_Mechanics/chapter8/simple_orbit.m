%simple_orbit.m - plots the orbit for a zero force case, u=1.5AU*sin(theta)=1/r
th=[0.35:0.01:pi-0.35];% angle range
r=1.5*csc(th);         % distance versus angle
subplot(1,2,1);plot(th,r);
axis([0 pi 0 5]),
xlabel('\theta (rad)','FontSize',14)
ylabel('r(\theta) (AU)','FontSize',14);
title('r(\theta) in the absence of a force','FontSize',14)
subplot(1,2,2);polar(th,r,'r')
title('Zero Force Straight Line Motion','FontSize',14)
xlabel('Distance in AU','FontSize',14)
