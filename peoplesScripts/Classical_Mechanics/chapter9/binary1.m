%binary1.m - binary star system given the eccentricity - using center of mass-
%relative coordinate concept 
clear;
m=5;m1=3;m2=m-m1;  %initial total, and individual masses (in Ms)
rcm=0; rmin=2.5;   %center of mass, rmin for relative coordinate in AU
e=0.6;             %given eccentricity (e=0 => circle)
%e=0.0;
th=[0:0.05:2*pi];
r=rmin*(1+e)./(1+e*cos(th));   %use orbit formula for reduced mass motion
rmax=rmin*(1+e)/(1-e);%use formula for rmax
%Note: reduced mass, m1 and m2, all have the same period - pick one
a=(rmin+rmax)/2;  %semimajor axis of the full reduced mass orbit
tau=sqrt(a^3/m);  %tau in units of tau_0 (tau_0=1 year)
str=cat(2,'Binary System, CM-Rel. Coord. Method: tau = ',num2str(tau,4),' tau_0,',...
' m_1 = ',num2str(m1,2),' m_s, m_2 = ',num2str(m2,2),' m_s');
r1=rcm-m2*r/m; r2=rcm+m1*r/m;%m1, m2 positions
[nr,nc]=size(r);             %get rows and columns of r
%polar plots are similar to reg plots of x, y with x=r.*cos(th); y=r.*sin(th), etc;
 for i=1:nc
     clf;
     hold on
     plot(0,rcm,'*')
     polar(th,r,'r:')        %show full paths
     polar(th,r1,'k--')
     polar(th,r2,'b-.')
     polar(th(i),r(i),'rs')  %place animation symbols
     polar(th(i),r2(i),'bo')
     polar(th(i),r1(i),'k.')
     polar(th(i),r2(i),'bo')
     pause(0.025)
 end
 xlabel('x (AU)','FontSize',14),ylabel('y (AU)','FontSize',14)
 title(str,'FontSize',11)
 h=legend('cm','\mu','m_1','m_2'); set(h,'FontSize',12)