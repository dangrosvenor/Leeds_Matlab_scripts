%conic3.m - script uses the conic section curve formula with rmim and
%eccentricity ye. Case of Rutherford scattering, given the initial 
%projectile energy and impact parameter.
clear; warning off;
rcm=0.0;                 %target position
za=2; zt=79;             %za=projectile, zt=target charges (2->alpha, 79->gold)
m=1;                     %projectile mass in units of alpha particle mass
vb=0.01965;              %velocity units (a_b/tau_b) (in inits of c=light speed)
ma=3730e6;               %alpha particle mass energy in eV
Ene=5e6;                 %initial projectile energy in eV
v0=sqrt(2*Ene/m/ma)/vb;  %initial speed in units of vb
b=20;                    %impact parameter
K=za*zt;                 %dimensionless force constant
ye=-sign(K)*sqrt(1+m^2*v0^4*b^2/K^2);% eccentricity
rmin=-m*v0^2*b^2/K/(1+ye);
L=100;%for window viewing purpose
% find theta limits for the asymptotes of r for this ye
% these occurs at the zeros of the denominator of the r(theta) equation
thmin=fzero(inline('ye*cos(x)+1'),[-pi/2,0],[],ye)+1.e-5;
thmax=fzero(inline('ye*cos(x)+1'),[0,pi/2],[],ye)-1.e-5;
THETA=pi-(thmax-thmin);%scattering angle
%----------------
%show the roots of the denominator 1+ye*cos(th) 
figure
x=[-2*pi:.1:2*pi];
plot(x,cos(x))
hold on
line([-2*pi,2*pi],[-1/ye,-1/ye],'Color','red')
plot(thmin,-1/ye,'ko',thmax,-1/ye,'ko')
text(thmin*(2+.3),(-1/ye)*(1+0.1),'\theta_{min}')
text(thmax*(1+0.1),(-1/ye)*(1+0.1),'\theta_{max}')
xlabel('\theta','FontSize',13),ylabel('cos\theta, -1/e','FontSize',13)
str=cat(2,' Asymptotes Angles (e = ',num2str(ye,3),'), \Theta=',num2str(THETA,3),' rad');
title(str,'FontSize',13)
h=legend('cos\theta','-1/e','Asymptotes',0);
set(h,'FontSize',13)
%----------------
th=[thmin:0.025:thmax];
r=rmin*(1+ye)./(1+ye*cos(th));
ths=th-thmin;% rotate orbital path counterclockwise by 'thmin'
             % to align asymptote with + x-axis                  
[nr,nc]=size(r);  %get rows and columns of r
figure
 for i=1:1:nc-10  %play all minus 10 to freeze particle
     clf;
     hold on
     plot(rcm,rcm,'*')
     polar(ths,r,'k:')        %show full path
     polar(ths(i),r(i),'b.')  %place animation symbols
     axis ([-L/2 L -L/2 L])
     line([-L/2,L],[0,0],'Color','red','LineStyle',':')
     pause(0.025)
 end
xlabel('x (a_b)','FontSize',13), ylabel('y (a_b)','FontSize',13)
str=cat(2,'r = r_{min}*(1+e)/(1+e*cos(\theta)) vs \theta; e = ',...
    num2str(ye,3),', r_{min}= ',num2str(rmin,3),' a_b, \Theta=',...
    num2str(THETA,3),' rad = ',num2str(THETA*360/2/pi,3),'^o');
str2=cat(2,'E_k= ',num2str(Ene,3),' eV',', b= ',num2str(b,3),' a_b',...
    ', m= ',num2str(m,3),' m_\alpha');
text(-L/2*(1-0.1),L*(1-0.1),str2)
title(str,'FontSize',11)
h=legend('Gold Target','Projectile Path','Projectile Particle',0);
set(h,'FontSize',13)
