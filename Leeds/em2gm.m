function [ phstoph lmstolm ] = em2gm(phis,lams,polphi,pollam)

% Constants
zrpi18  = 57.2957795;
zpir18  = 0.0174532925;

% General settings
zsinpol = sin(zpir18*polphi);
zcospol = cos(zpir18*polphi);
zlampol = zpir18*pollam;
zphis    = zpir18*phis;
zlams    = lams;
if (zlams > 180) 
   zlams = zlams - 360;
end
zlams    = zpir18*zlams;

% Transform latitude
zarg    = zcospol*cos(zphis).*cos(zlams) + zsinpol*sin(zphis);
phstoph = zrpi18*asin(zarg);

% Transform longitude
zarg1   = sin(zlampol)*(-zsinpol*cos(zlams).*cos(zphis)  + ...
                         zcospol*            sin(zphis)) - ...
          cos(zlampol)*          sin(zlams).*cos(zphis);
zarg2   = cos(zlampol)*(-zsinpol*cos(zlams).*cos(zphis)  + ...
                         zcospol*            sin(zphis)) + ...
          sin(zlampol)*          sin(zlams).*cos(zphis);
lmstolm = zrpi18*atan2(zarg1, zarg2);
