%curl_ex.m
%This script plots and evaluates the function V=-y i + x j 
%and its curl. Other functions are possible as well.
warning off; %supress unwanted warnings by plotter, if needed
clear; vmin=-1; vmax=1;
xmax=vmax; ymax=vmax; zmax=vmax; % x,y,z upper limits
xmin=vmin; ymin=vmin; zmin=vmin;
vs=0.3; xs=vs; ys=vs; zs=vs; % step size
N=(vmax-vmin)/vs; % number of points to be plotted is NxNxN
m=round(N/2+1); zm=zmin+(m-1)*zs; % z at which the vector function is plotted
%[x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);
[x,y,z]=meshgrid(xmin:xs:xmax,ymin:ys:ymax,zmin:zs:zmax);
% The desired vector is F = fx i + fy j + fz k. Uncomment appropriate example
% Example 1
  fx=-y;fy=x;fz=0.*z;str1='Vector (-y,x,0)';
% Example 2
% fx=x.*z;fy=-y.^2;fz=2*x.^2.*y;str1='Vector (x*z,-y^2,2x^2*y)';
% Example 3
% fx=y.^2;fy=(2*x.*y+z.^2);fz=2*y.*z;str1='Vector (y^2,(2*x*y+z^2),2*y*z)';
% Example 4
% fx=x.*y;fy=y.*z;fz=z.*x;str1='Vector (xy,yz,zx)';
% Example 5
% fx=x.^2.*y;fy=y.^2.*z;fz=z.^2.*x;str1='Vector (x^2*y,y^2*z,z^2*x)';
% quiver3(x,y,z,fx,fy,fz,2); %uses scaling
 quiver3(x,y,z,fx,fy,fz); %draws arrows in three dimensions
 title (str1,'FontSize',14)
 xlabel('x','FontSize',14), ylabel('y','FontSize',14)
 zlabel('z','FontSize',14)
 [curlx,curly,curlz,cav] = curl(x,y,z,fx,fy,fz);%Curl, angular vel. in rad/sec of F
figure
% surf(x(:,:,m),y(:,:,m),cav(:,:,m)) % surface plot of the angular velocity
%================================
 x2=x(:,:,m); y2=y(:,:,m); fx2=fx(:,:,m); fy2=fy(:,:,m);
 cav2 = curl(x2,y2,fx2,fy2);% curl's average angular velocity in one plane of the volume
 surf(x2,y2,cav2) % surface plot of the angular velocity
 shading interp
 hold on;
 quiver(x2,y2,fx2,fy2);%plots the function vectors at z=zm
 str2=cat(2,'\omega_z surface and function vectors at ','z=',num2str(zm,3));
 str4=cat(2,'\omega_z, f_x, f_y, ');
 title (str2,'FontSize',14)
 xlabel('x','FontSize',14), ylabel('y','FontSize',14)
 zlabel(str4,'FontSize',14)
 hold off
%================================
figure
quiver3(x,y,z,curlx,curly,curlz); %draws curl's result as arrows in three dimensions
str3=cat(2,'Curl');
title (str3,'FontSize',14) 
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
zlabel('z','FontSize',14)
% %subplot(2,2,4) % Lines below prove equality to cav2 = curl(x2,y2,fx2,fy2)/2
% figure
% cav3=0.5*curlz(:,:,m);% same as curl's average angular velocity about the axis
% surf(x2,y2,cav3);%plots the velocity vectors
% shading interp
% hold on;
% quiver(x2,y2,fx2,fy2);%plots the velocity vectors
% hold off
