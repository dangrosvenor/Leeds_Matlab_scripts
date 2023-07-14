%ellipse.m - script to draw an ellipse with minimum radius rmim and
%eccentricity e.
clear;
rmin=2.5;
e=0.6;
th=[0:0.1:2*pi];
r=rmin*(1+e)./(1+e*cos(th));
rmax=rmin*(1+e)./(1+e*cos(pi));
polar(th,r)
str=cat(2,'Ellipse with (rmin, rmax) = (',num2str(rmin,3),',',...
    num2str(rmax,3),')AU, e = ',num2str(e,3));
title(str,'FontSize',14)