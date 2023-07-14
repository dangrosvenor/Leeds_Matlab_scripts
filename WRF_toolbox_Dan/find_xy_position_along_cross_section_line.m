function [x,y]=find_xy_position_along_cross_section_line(x_line,y_line,d)
%finds the the x,y position for the point d km along the cross section line
%described by x_line and y_line

c = diff(y_line)/diff(x_line);
dy = sign(diff(y_line)) * sqrt(d^2/(1+1/c^2));
dx = dy/c;

x=x_line(1)+dx;
y=y_line(1)+dy;
