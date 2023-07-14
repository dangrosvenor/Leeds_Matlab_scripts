%gradient_ex.m
%This script plots and evaluates the function f=x*exp(-(x^2+y^2+z^2)) and its gradient
%The plots are done versus x,y at a certain z value
warning off; %supress unwanted warnings by plotter if needed
clear;
vmax=2.0;
xmax=vmax; ymax=vmax; zmax=vmax; % x,y,z limits
vs=0.1;
xs=vs; ys=vs; zs=vs; % step size
N=2*vmax/vs; % number of points to be plotted is NxNxN
dv=0.1;
dx=dv; dy=dv; dz=dv; % used in the gradient
m=round(N/2+5); zm=-zmax+(m-1)*zs; % value of z at which we plot f(x,y,z)
[x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);
f=x.*exp(-(x.^2+y.^2+z.^2)); % the desired function
[dfx,dfy,dfz] = gradient(f,dx,dy,dz);%gradient of f(x,y,z)
%mesh(x(:,:,m),y(:,:,m),f(:,:,m)) % can do a mesh if desired
surf(x(:,:,m),y(:,:,m),f(:,:,m)) % surface plot
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
zlabel('f(x,y,z)','FontSize',14)
str=cat(2,'f(x,y,z)=x*exp(-x^2-y^2-z^2) at ','z=',num2str(zm,3));
title(str,'FontSize',14)
figure
%contour(x(:,:,m),y(:,:,m),f(:,:,m),20)%contour plot 20 line case
[C,h] = contour(x(:,:,m),y(:,:,m),f(:,:,m));% generate contour plot
clabel(C,h,'FontSize',12)% add countour labels
hold on
quiver(x(:,:,m),y(:,:,m),dfx(:,:,m),dfy(:,:,m))%draw the gradient as arrows at x
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
str=cat(2,'f(x,y,z)=x*exp(-x^2-y^2-z^2) at ','z=',num2str(zm,3));
title(str,'FontSize',14)