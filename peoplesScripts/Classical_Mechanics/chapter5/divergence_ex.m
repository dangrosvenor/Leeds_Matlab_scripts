%divergence_ex.m
%This script plots and evaluates the function r=x i + y j + z k 
%and its divergence. Other functions are possible as well.
warning off; %supress unwanted warnings by plotter, if needed
clear; vmin=-1; vmax=1;
xmax=vmax; ymax=vmax; zmax=vmax; % x,y,z upper limits
xmin=vmin; ymin=vmin; zmin=vmin;
vs=0.4; xs=vs; ys=vs; zs=vs; % step size
N=(vmax-vmin)/vs; % number of points to be plotted is NxNxN
m=round(N); zm=zmin+(m-1)*zs; % z at which the vector function is plotted
%[x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);
[x,y,z]=meshgrid(xmin:xs:xmax,ymin:ys:ymax,zmin:zs:zmax);
% The desired vector is F = fx i + fy j + fz k. Uncomment appropriate example
% Example 1
fx=x; fy=y; fz=z;str1='Vector r';
% Example 2
%fx=x.*z; fy=-y.^2; fz=2*x.^2.*y;str1='Vector (x*z,-y^2,2*x^2*y)';
% Example 3
% fx=x.^3; fy=y; fz=z;str1='Vector (x^3,y,z)';
% Example 4
% fx=x.^2.*y; fy=y.^2.*z; fz=z.^2.*x;str1='Vector (x^2*y,y^2*z,z^2*x)';
%quiver3(x,y,z,fx,fy,fz,2); %uses scaling
quiver3(x,y,z,fx,fy,fz); %draws arrows in three dimensions
title (str1,'FontSize',14); 
xlabel('x','FontSize',14); ylabel('y','FontSize',14); 
zlabel('z','FontSize',14)
div = divergence(x,y,z,fx,fy,fz);%divergence of F
figure
surf(x(:,:,m),y(:,:,m),div(:,:,m)) % surface plot
str2=cat(2,'Divergence at ','z=',num2str(zm,3));
title(str2,'FontSize',14) 
xlabel('x','FontSize',14), ylabel('y','FontSize',14)