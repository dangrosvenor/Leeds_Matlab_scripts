function [xx,yy]=rotate(x,y,thi,Y,xmin,ymin,ymax)

yd=Y*(y-ymin)/(ymax-ymin);
xd=x-xmin; %transform x and y so that are in "page" co-ordinates, i.e. y scale is simlar to x scale
%h=x.*sin(thi);
%xx=xd.*cos(thi);
%yy=yd-h;

%xx=xd;

%yy=-xd.*sin(thi).*(1-cos(thi));
% b=yd-(xd./cos(thi));
% yy=(b./sin(thi)) + xd.*tan(thi);
% rr=( ((xd).^2)+((yd).^2) ).^0.5;
% angle=atan( (yd)./(xd) );
% xx=rr.*cos(angle-thi);
% yy=(rr.*sin(angle-thi));

yy=yd.*sin(thi);
xx=xd+(yd.*cos(thi));