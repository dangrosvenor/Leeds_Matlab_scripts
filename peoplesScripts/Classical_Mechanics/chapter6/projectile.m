%projectile.m - without air resistance
%plots the projectile trajectory z(x)=viz*x/vix-0.5*g*(x/v0x)^2
clear;
vi=5; zi=0; xi=0; g=9.8;%initial: speed, positions, gravity
xmin=xi;xs=0.005; thi=15; ths=15; thm=75;% x, angle ranges
hold on;
for th=thi:ths:thm;
    vix=vi*cos(th*2*pi/360);% initial speed components
    viz=vi*sin(th*2*pi/360);
    xmax=2*vix*viz/g+xi;
    N=(xmax-xmin)/xs+1; M=round(N/2);
    x=[xmin:xs:xmax];
    z=zi+viz*(x-xi)/vix-0.5*g*((x-xi)/vix).^2;
    plot(x,z);
str1=cat(2,'\theta =',num2str(th),'^o');
text(x(M+25),z(M+25)+0.025,str1,'FontSize',10,'Color','red');
end
str3=cat(2,'Trajectory plots for different initial angle');
title(str3,'FontSize',14); ylabel('z(x)','FontSize',14)
xlabel('x','FontSize',14); grid on;
