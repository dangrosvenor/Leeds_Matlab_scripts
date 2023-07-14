function [u,v]=uv_from_winddir(speed,winddir)
%get u and v components from teh wind direction
%[u,v]=uv_from_winddir(speed,winddir)
% speed in m/s, windir in degrees from north (clockwise)

u = speed.*sin(winddir.*pi/180);
v = speed.*cos(winddir.*pi/180);

