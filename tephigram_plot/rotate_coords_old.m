function [x_new,y_new]=rotate_coords(x,y,theta);
% rotates
%[th,r]=cart2pol(x,y);
%th=th-theta;
%[x_new,y_new]=pol2cart(th,r);

%x_new=cos(-theta).*x -sin(-theta).*y;
%y_new=sin(-theta).*x+ cos(-theta).*y;
x_new=x;
y_new=y;

% nb only want to rotate center point?
