%v_and_f.m
% plots V=(1/x^3-1/x^2) and the related force F=3/x^4-2/x^3
clear;
xb=3/2;              %the bond length
vmin=-4/27;          %potential min determines size of plotting window
xmin=0.5;            %minimum value of x
xmax=5*xb;           %use the bond length to determine the maximum x 
NPTS=100;            %number of points
x=[xmin:2/NPTS:2*xmax]; %distance array
V=1./x.^3-1./x.^2;   %the potential array
F=3./x.^4-2./x.^3;   %the force array
plot(x,V,'b:',x,F,'r-.');
line([0,xmax],[0,0],'Color','k','LineStyle','--') %draw a zero line
axis([xmin, xmax,1.5*vmin,-1.5*vmin])
title('Plot of V(x)=(1/x^3-1/x^2) and F=3/x^4-2/x^3 vs position','FontSize',14)
ylabel('V(x), F(x) (u_0)','FontSize',14);
xlabel('x(a_0)','FontSize',14);
h=legend('V','F',0); set(h,'FontSize',14)
