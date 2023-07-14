function [xx,yy]=rotate(x,y,thi)

rr=( ((x).^2)+((y).^2) ).^0.5;
angle=atan( (y)./(x) );
dan=size(thi)
xx=rr.*cos(angle-thi);
yy=(rr.*sin(angle-thi));