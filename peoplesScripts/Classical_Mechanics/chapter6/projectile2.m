%projectile2.m - motion with air resistance
%plots the projectile trajectory:
%z(x)=(m*g/c+viz)*(x-xi)/vix+m^2*g*log(1-c*(x-xi)/m/vix)/c^2
clear;
vi=25;zi=0;xi=0;th=45;g=9.8;%initial: speed, positions, angle, gravity
c=0.5; m=1; % Drag Coeff and mass
vix=vi*cos(th*2*pi/360); viz=vi*sin(th*2*pi/360);%initial vi components
xmin=xi;xs=0.005;xmax=2*vix*viz/g;% Range without air resistance
x=[xmin:xs:xmax];
z=zi+viz*(x-xi)/vix-0.5*g*((x-xi)/vix).^2;% z(x) without air resistance
tmax=2*viz/g; %no air resistance round trip time for initial guess
%fprintf('Free Motion: t_R %f sec, R= %f m\n',tmax,xmax)
stra=cat(2,'t_R=',num2str(tmax,4),',R=',num2str(xmax,4));
%======= Newton Raphson Method for time to land
Nt=0; tdiff=1;
while Nt < 5 & tdiff > 0.0001, % iteration conditions
    Nt=Nt+1;                   % counter
    f=(m/c+viz/g)*(1-exp(-c*tmax/m))-tmax;% function whose zero we seek
    fp=c*(m/c+viz/g)*exp(-c*tmax/m)/m-1;  % derivative 
    tdiff=f/fp;tmax=tmax-tdiff;% tmax becomes the new guess
end
% ============================
xmaxf=m*vix*(1-exp(-c*tmax/m))/c;%Range with air resistance
xf=[xmin:xs:xmaxf];
zf=(m*g/c+viz)*(xf-xi)/vix+m^2*g*log(1-c*(xf-xi)/m/vix)/c^2;
plot(x,z,'b--',xf,zf,'r');
str0=cat(2,'\theta =',num2str(th),'^o',',c=',...
    num2str(c),'kg/s', ',vi=',num2str(vi),'m/s');
str1=cat(2,'Trajectory plots w/o Air Resistance',' at ',str0);
title(str1,'FontSize',13); ylabel('z(x)','FontSize',14)
xlabel('x','FontSize',14); grid on;
h=legend('Free Motion','Air Resistance',1); set(h,'FontSize',12)
%fprintf('With Air Resistance: t_R %f sec, R= %f m\n',tmax,xmaxf)
strb=cat(2,'t_R=',num2str(tmax,4),',R=',num2str(xmaxf,4));
Ma=round(4*length(x)/5); Mb=round(4*length(xf)/5);
text(x(Ma)*(1+0.01),z(Ma),stra,'FontSize',9,'Color','blue');
text(xf(Mb)*(1+0.02),zf(Mb),strb,'FontSize',9,'Color','red');
