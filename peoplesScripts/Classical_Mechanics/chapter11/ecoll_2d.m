%ecoll_2d.m - two-dimensional collision script to find the final velocities
%of two particles and their scattering angles.
%needed input is the initial velocities, their masses, and one of the 
%scattered angles. This version was made from 'ecoll_2d0.m', but it's based
%on the use of the mmfsolve.m if the optimization toolbox in not available
clear;
global v1i v2i m1 m2 tp10 tp20 tp1
%th10=0.0;th20=0.0;th1=35;m1=0.5; m2=0.75; v1i=2; v2i=0.0;%example of m2 at rest
%th10=20;th20=10;th1=15;m1=0.5; m2=0.75; v1i=2; v2i=1.5;%masses, initial angles,speeds
th10=40;th20=20;th1=30;m1=0.5; m2=1.5; v1i=1.5; v2i=0.75;%masses, init angles,speeds
tp10=th10*pi/180;tp20=th20*pi/180;tp1=th1*pi/180;      %convert angles to radians
v1fg=v1i/2;v2fg=v1i*sqrt(0.75*m1/m2);tp2g=pi/2-tp1;    %guesses - change as needed
%opts =mmfsolve('Display','off');                      %Turn off Display
opts =mmfsolve('FunTol',1e-7,'MaxIter',100);           %alternate options
X0=[v1fg;v2fg;tp2g;];
X=mmfsolve(@ecoll_2dfun,X0,opts);
v1f=X(1);v2f=X(2);tp2=X(3);th2=tp2*360/2/pi;
%------------------- plotting ------------------------------------------
line([-v1i*cos(tp10),0],[v1i*sin(tp10),0],'Color','k')                 %initial p m1
line([-v2i*cos(tp20),0],[-v2i*sin(tp20),0],'Color','m','LineStyle',':')%initial p m2
line([0,v1f*cos(tp1)],[0,v1f*sin(tp1)],'Color','b','LineStyle','--')   %scattered m1
line([0,v2f*cos(tp2)],[0,-v2f*sin(tp2)],'Color','r','LineStyle','-.')  %scattered m2
h=legend('m_1 initial','m_2 initial','m_1 final','m_2 final',4);
set(h,'FontSize',11)
va=max(v1i,v2i);vb=max(v1f,v2f);                       %pick higher values for axis
line([-va,vb],[0,0],'LineStyle',':')                   %shows x axis
line([0,0],[-va,va],'LineStyle',':')                   %shows y axis
axis([-va va -vb vb])
xlabel('p_x','FontSize',14)
ylabel('p_y','FontSize',14)
str0=cat(2,'Two mass elastic collision: m1, m2 = ',num2str(m1,3),...
           ', ',num2str(m2,3),' kg');
str1=cat(2,'v_{1i}=',num2str(v1i,3),'m/s, v_{2i}=',num2str(v2i,3),...
        'm/s, \theta_{10}=',num2str(th10,3),'^o, \theta_{20}=',...
        num2str(th20,3),'^o');
str2=cat(2,'v_{1f}=',num2str(v1f,3),'m/s, v_{2f}=',num2str(v2f,3),...
        'm/s, \theta_1=',num2str(th1,3),'^o, \theta_2=', num2str(th2,3),'^o');
title(str0,'FontSize',14);
text(-va*(1-.05),vb*(1-0.05),str1,'FontSize',[11]);
text(-va*(1-.05),vb*(1-0.175),str2,'FontSize',[11]);
text(-v1i*cos(tp10)*(1-0.5),v1i*sin(tp10)*(1-.85),'\theta_{10}','FontSize',[9])
if v2i~= 0,
text(-v2i*cos(tp20)*(1-0.5),-v2i*sin(tp20)*(1-0.85),'\theta_{20}','FontSize',[9])
end
text(v1f*cos(tp1)*(1-0.75),v1f*sin(tp1)*(1-0.92),'\theta_1','FontSize',[9])
text(v2f*cos(tp2)*(1-0.75),-v2f*sin(tp2)*(1-0.92),'\theta_2','FontSize',[9])
%------------------- test results, Momentum, Energy -----------------------------
Ei=(m1*v1i^2+m2*v2i^2)/2;
Ef=(m1*v1f^2+m2*v2f^2)/2;
Pxi=m1*v1i*cos(tp10)+m2*v2i*cos(tp20);
Pxf=m1*v1f*cos(tp1)+m2*v2f*cos(tp2);
Pyi=-m1*v1i*sin(tp10)+m2*v2i*sin(tp20);
Pyf=m1*v1f*sin(tp1)-m2*v2f*sin(tp2);
fprintf('Momentum, Energy Check: m1=%7.4f, m2=%7.4f\n',m1,m2)
fprintf('v1i=%7.4f, v2i=%7.4f, th10=%7.4f, th20=%7.4f\n',v1i,v2i,th10,th20)
fprintf('v1f=%7.4f, v2f=%7.4f, th1=%7.4f, th2=%7.4f\n',v1f,v2f,th1,th2)
fprintf('Ei=%7.4f, Ef=%7.4f, Pxi=%7.4f Pxf=%7.4f Pyi=%7.4f Pyf=%7.4f\n',...
         Ei,Ef,Pxi,Pxf,Pyi,Pyf);
