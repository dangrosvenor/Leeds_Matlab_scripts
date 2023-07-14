function [dlon]=distlatlon_inverse(d,latA,latB)
%function [d,aob]=distlatlon(latA,lonA,latB,lonB)
%d is distance in km
%also returns aob in degrees
%aob is therfore the angle subtended by the great circle (i.e. the
%angle between the two points with the Earth's centre as the origin
latA=latA *pi/180;
latB=latB *pi/180;
%lonA=lonA *pi/180;
%lonB=lonB *pi/180;

R=6378.140; %radius of earth (km)

% aob=acos( cos(latA).*cos(latB).*cos(lonB-lonA) + sin(latA).*sin(latB) );
% %here aob is therfore the angle subtended by the great circle (i.e. the
% %angle between the two points with the Earth's centre as the origin
% d=R*aob;
% aob=aob*180/pi;

%Invert the above to get lonB-lonA
aob = d./R;
dlon = acos( (cos(aob) - sin(latA)*sin(latB) ) / (cos(latA)*cos(latB)) );
dlon = dlon*180/pi;






