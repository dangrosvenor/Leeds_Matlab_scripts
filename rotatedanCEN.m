function [xx,yy]=rotate(x,y,thi)

rr=( ((x-230).^2)+((y-2.1).^2) ).^0.5;
angle=atan( (y-2.1)./(x-230) );
xx=rr.*cos(angle-thi);
yy=(rr.*sin(angle-thi));