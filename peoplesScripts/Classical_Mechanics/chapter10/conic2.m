%conic2.m - script for conic section projectile hyperbolic orbit
%with negative eccentricity ye.
clear; warning off;
ye=-1.2; %eccentricity
rmin=1;
L=5;%for window viewing purpose
% find theta limits for the assymptotes of r for this ye
% these occur at the zeros of the denominator of the r(theta) equation
% opt=optimset('TolX',1.0e-10); %optional convergence criteria
% thmin=fzero(inline('ye*cos(x)+1'),[-pi/2,0],opt,ye)+1.e-3;
% thmax=fzero(inline('ye*cos(x)+1'),[0,pi/2],opt,ye)-1.e-3;
thmin=fzero(inline('ye*cos(x)+1'),[-pi/2,0],[],ye)+1.e-5;
thmax=fzero(inline('ye*cos(x)+1'),[0,pi/2],[],ye)-1.e-5;
THETA=pi-(thmax-thmin);%scattering angle
%--------show asymptote angles----------
figure%plots the cos(x) versus -1/ye to show assymptote angles
x=[-2*pi:.1:2*pi];
plot(x,cos(x),'b--')
hold on
line([-2*pi,2*pi],[-1/ye,-1/ye],'Color','red')
plot(thmin,-1/ye,'ko',thmax,-1/ye,'ko')
text(thmin*(2+.3),(-1/ye)*(1+0.1),'\theta_{min}')
text(thmax*(1+0.1),(-1/ye)*(1+0.1),'\theta_{max}')
xlabel('\theta','FontSize',14),ylabel('cos\theta, -1/e','FontSize',14)
str=cat(2,' Asymptotes Angles (e = ',num2str(ye),'), \Theta=',num2str(THETA),' rad');
title(str,'FontSize',12)
h=legend('cos\theta','-1/e','Asymptotes',4); set(h,'FontSize',12)
%----------------------------------------
th=[thmin:0.001:thmax]; %range of orbit plot
r=rmin*(1+ye)./(1+ye*cos(th));
[nr,nc]=size(r);             %get rows and columns of r
%polar plots are similar to reg plots of x, y with x=r.*cos(th); y=r.*sin(th), etc;
figure
for i=1:25:nc
     clf;
     hold on
     plot(0,0,'*')
     polar(th,r,'r:')        %show full path
     polar(th(i),r(i),'r.')  %place animation symbols
     axis ([-L/L L -L L])
     pause(0.025)
 end
xlabel('x (Arbitrary Units)','FontSize',14)
ylabel('y (Arbitrary Units)','FontSize',14)
str=cat(2,'Orbit: r = r_{min}*(1+e)/(1+e*cos(\theta) vs \theta for e = ',...
    num2str(ye),', r_{min}=',num2str(rmin),', \Theta=',num2str(THETA),' rad');
title(str,'FontSize',12)
h=legend('Target','Projectile'); set(h,'FontSize',12)
 
