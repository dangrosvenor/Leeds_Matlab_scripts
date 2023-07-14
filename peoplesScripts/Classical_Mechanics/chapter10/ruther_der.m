%ruther_der.m: returns the derivatives for Rutherford scattering trajectory
function derivs = ruther_der( t, w, flag,K,m)
% za=projectile charge, zt=charge of target, m=projectile mass
% Entries in the vector of dependent variables are:
% w(1,2,...8)=x1,v1x,y1,v1y,x2,v2x,y2,v2y
wr=sqrt(w(1).^2 + w(3).^2);
derivs = [w(2); K*w(1)./wr.^3/m; w(4);K*w(3)./wr.^3/m];
