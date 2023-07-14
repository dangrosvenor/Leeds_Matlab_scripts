%ecoll_2d0.m - two dimensional collision script to find the final velocities
%of two particles and their scattering angles.
%needed input is the initial velocities, their masses, and one of the 
%scattered angles.
%Based on the use of the fsolve.m function in the optimization toolbox
clear;
global v1i v2i m1 m2 tp10 tp20 tp1
%th10=0.0;th20=0.0;th1=35;m1=0.5; m2=0.75; v1i=2; v2i=0.0;%example of m2 at rest
%th10=20;th20=10;th1=15;m1=0.5; m2=0.75; v1i=2; v2i=1.5;%masses, initial angles,speeds
th10=40;th20=20;th1=30;m1=0.5; m2=1.5; v1i=1.5; v2i=0.75;%masses, initial angles,speeds
tp10=th10*pi/180;tp20=th20*pi/180;tp1=th1*pi/180;      %convert angles to radians
v1fg=v1i/2;v2fg=v1i*sqrt(0.75*m1/m2);tp2g=pi/2-tp1;    %guesses - change as needed
X0=[v1fg,v2fg,tp2g];
opts = optimset('Display','off');                              %Turn off Display
%opts=optimset('MaxIter',2000,'MaxFunEvals',20000,'TolX',1e-8); %alternate options
%X = fsolve(@ecoll_2dfun0,X0,opts,v1i,v2i,m1,m2,tp10,tp20,tp1);%this without a "global" statement
X = fsolve(@ecoll_2dfun0,X0,opts);                             %this with a "global" statement
v1f=X(1);v2f=X(2);tp2=X(3);th2=tp2*360/2/pi;
%------------------- plotting ------------------------------------------
line([-v1i*cos(tp10),0],[v1i*sin(tp10),0],'Color','k') %initial momentum particle 1
line([-v2i*cos(tp20),0],[-v2i*sin(tp20),0],'Color','m')%initial momentum particle 2
line([0,v1f*cos(tp1)],[0,v1f*sin(tp1)],'Color','b')    %scattered particle 1
line([0,v2f*cos(tp2)],[0,-v2f*sin(tp2)],'Color','r')   %scattered particle 2
legend('m_1 initial','m_2 initial','m_1 final','m_2 final',4)
va=max(v1i,v2i);vb=max(v1f,v2f);                       %pick higher values for axis
line([-va,vb],[0,0],'LineStyle',':')                   %shows x axis
line([0,0],[-va,va],'LineStyle',':')                   %shows y axis
axis([-va va -vb vb])
xlabel('p_x')
ylabel('p_y')
str0=cat(2,'Two mass elastic collision: m1, m2 = ',num2str(m1,3),...
           ', ',num2str(m2,3),' kg');
str1=cat(2,'v_{1i}=',num2str(v1i,3),'m/s, v_{2i}=',num2str(v2i,3),...
        'm/s, \theta_{10}=',num2str(th10,3),'^o, \theta_{20}=',...
        num2str(th20,3),'^o');
str2=cat(2,'v_{1f}=',num2str(v1f,3),'m/s, v_{2f}=',num2str(v2f,3),...
        'm/s, \theta_1=',num2str(th1,3),'^o, \theta_2=', num2str(th2,3),'^o');
title(str0);
text(-va*(1-.05),vb*(1-0.05),str1,'FontSize',[10]);
text(-va*(1-.05),vb*(1-0.175),str2,'FontSize',[10]);
text(-v1i*cos(tp10)*(1-0.5),v1i*sin(tp10)*(1-.85),'\theta_{10}','FontSize',[8])
if v2i~= 0,
text(-v2i*cos(tp20)*(1-0.5),-v2i*sin(tp20)*(1-0.85),'\theta_{20}','FontSize',[8])
end
text(v1f*cos(tp1)*(1-0.75),v1f*sin(tp1)*(1-0.92),'\theta_1','FontSize',[8])
text(v2f*cos(tp2)*(1-0.75),-v2f*sin(tp2)*(1-0.92),'\theta_2','FontSize',[8])
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
