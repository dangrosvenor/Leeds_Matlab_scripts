%pend0.m
%This program shows the relationship between A1 and the initial angle, for the
%non-linear approximation of the pendulum's sin(theta) term
clear;
w0=1;
a3=w0^2/6;
imax=10;
tol=1.e-5;
th0=0;
thmax=90;
N=25;
dth=(thmax-th0)/N;
for j=1:N,
  th(j)=th0+(j-1)*dth;
  x=th(j)-1;%initial guess
  xn=999;
  f=999;
  i = 0;
  while (abs(xn-x) >= tol) & (f ~= 0.0) & (i < imax)
      x=xn;
      f=x+a3*x^3/(27*a3*x^2-32*w0^2)-th(j);
      fp=1+3*a3*x^2/(27*a3*x^2-32*w0^2)...
          -54*a3^2*x^4/((27*a3*x^2-32*w0^2)^2);
      xn=x-f/fp;
      i = i + 1;
  end
A1(j)=xn;
end
plot(th,A1)
xlabel('\theta_0 (Degrees)','FontSize',14)
ylabel('A_1 (Degrees)','FontSize',14)
title('Relationship between A_1 and \theta_0','FontSize',14)
