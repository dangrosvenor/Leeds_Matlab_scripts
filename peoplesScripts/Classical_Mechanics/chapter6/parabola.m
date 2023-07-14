%parabola.m
%plots parabolas (x-x0)^2=4e(z-z0), whose curvatures=1/2e,
%(x0,z0) is the vertex, and (x0,z0+e) is the focal point position
clear;
x0=1; z0=1; %parameters
xmin=-3; xmax=3; xs=0.1; % x-range
x=[xmin+x0:xs:xmax+x0];  % x array
hold on; zt=0; zb=0; % useful in plotting
for e=-0.5:-0.5:-1.5, %vary e
z=z0+(x-x0).^2/4/e;
plot(x,z,x0,z0,'r.',x0,z0+e,'rx');
str1=cat(2,'F (',num2str(e),')');% focal point
str2=cat(2,'e=',num2str(e));
text(x0*(1+0.2),z0*(1+0.1)+e,str1,'FontSize',12,'Color','red');
text(x(4),z(4),str2,'FontSize',10,'Color','red');
if abs(max(z)) > abs(zt); zt=max(z); end
if abs(min(z)) > abs(zb); zb=min(z); end
end
str3=cat(2,'(x-x_0)^2=4e(z-z_0)',' with x_0=',num2str(x0),...
   ',z_0=',num2str(z0));
%title('Plot of (x-x_0)^2=4e(z-z_0)')
title(str3,'FontSize',10); ylabel('z(x)','FontSize',14)
xlabel('x','FontSize',14); grid on;
axis([xmin+x0 xmax+x0 zb*(1-0.2) zt*(1+0.2)]);
text(x0*(1+0.2),z0*(1+0.1),'Vertex','FontSize',12,'Color','red');
