function [xd,yd]=pagetrans(x,y,Y,xmin,ymin,ymax)

yd=Y*(y-ymin)/(ymax-ymin);
xd=x-xmin;