function [xx,yy]=rotate(x,y,thi,Y,xmin,ymin,ymax)

yd=Y*(y-ymin)/(ymax-ymin);
xd=x-xmin; %transform x and y so that are in "page" co-ordinates, i.e. y scale is simlar to x scale
rr=( ((xd).^2)+((yd).^2) ).^0.5;
angle=atan( (yd)./(xd) );
xx=rr.*cos(angle-thi);
yy=(rr.*sin(angle-thi));