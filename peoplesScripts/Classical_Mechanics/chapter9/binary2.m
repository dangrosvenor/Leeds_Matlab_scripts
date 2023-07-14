%binary2.m - program to simulate the binary mass system in a fully numerical way
clear;
rcm=0.0;             %center of mass
c=4*pi^2;            %constant
m=5;m1=3;m2=m-m1;    %initial total, and individual masses (in Ms)
tmax=6.988;          %maximum time in units of tau_0 (tau_0=1 yr)
%tmax=1.77;
x10=-1;y10=0;x20=-m1*x10/m2;y20=0;
x00=x20-x10;         %relative coordinate circular orbit parameter
e=0.6;               %given eccentricity (e=0 => circle)
%e=0.0;   
d=sqrt(1+e); %e > 0 =>deviation from circular orbit
v1x0=0;v1y0=d*(-m2*2*pi*sqrt(m/x00)/m);v2x0=0;v2y0=-m1*v1y0/m2;
ic1=[x10,v1x0,y10,v1y0,x20,v2x0,y20,v2y0]; %initial conditions
%Use MATLAB's Runge-Kutta (4,5) formula (uncomment/comment as needed)
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-4);%user set Tolerances
[t,w]=ode45('binary2_der',[0.0,tmax],ic1,opt,c,m1,m2);%with set tolerance
%[t,w]=ode45('binary2_der',[0.0,tmax],ic1,[],c,m1,m2);%with default tolerance
%calculate period obtained in units of tau_0
r1=sqrt(w(:,1).^2+w(:,3).^2); r2=sqrt(w(:,5).^2+w(:,7).^2);%get radii
rmin1=min(r1); rmax1=max(r1); rmin2=min(r2); rmax2=max(r2);%get min, max r's
a1=(rmin1+rmax1)/2; a2=(rmin2+rmax2)/2; a=a1+a2;%get the a's
tau=sqrt(a^3/m); %get the period
str=cat(2,'Binary Syst. Numeric Method: tau = ',num2str(tau,4),' tau_0, ',...
' m_1 = ',num2str(m1,2),' m_s, m_2 = ',num2str(m2,2),' m_s');
[nr,nc]=size(w);      %get rows and columns of r
figure
for i=1:nr
    clf;
    plot(rcm,rcm,'*')% center of mass
    hold on
    plot(w(:,1),w(:,3),'k:',w(:,5),w(:,7),'b-.')
    plot(w(i,1),w(i,3),'k.','MarkerSize',20)
    plot(w(i,5),w(i,7),'bo','MarkerSize',8)
    axis([-6.0 4.0 -3.0 3.0]);
    pause(0.0125)
end
xlabel('x (AU)','FontSize',14),ylabel('y (AU)','FontSize',14)
title(str,'FontSize',12)
h=legend('cm','m_1','m_2'); set(h,'FontSize',12)
figure
subplot(1,2,1)%plot radii vs time
plot(t,r1,'k.',t,r2,'bo')
xlabel('t (yr)','FontSize',14),ylabel('r_1,r_2 (AU)','FontSize',14)
title('Binary Syst. Numeric: r versus t ','FontSize',12)
subplot(1,2,2)%plot speeds vs time
v1=sqrt(w(:,2).^2+w(:,4).^2); v2=sqrt(w(:,6).^2+w(:,8).^2);%get speeds
plot(t,v1,'k.',t,v2,'bo')
xlabel('t (yr)','FontSize',14),ylabel('v_1,v_2 (AU/yr)','FontSize',14)
title('Binary Syst. Numeric: v versus t ','FontSize',12)
