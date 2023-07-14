function [range, d] = dist_noctilucent(angle_elev_deg,H)

angle_elev = pi/180 * angle_elev_deg; %convert to radians

R=6371; %km

A=R+H;

%angle between observer and nadir for the clouds from sin formula
r = asin( R.* sin( pi/2 + angle_elev ) ./ A );

theta = pi/2 - angle_elev - r;
range = R.*theta; %arc distance long earth's surface to below cloud.
d = sqrt( R.^2 + A.^2 - 2.*A.*R.*cos(theta) );
