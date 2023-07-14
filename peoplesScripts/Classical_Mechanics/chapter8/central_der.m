%central_der.m: returns the derivatives for -a*r^p central force
function derivs = central_der( t, w, flag,a,L,m,p)
% a=force strength, L=angular momentum, m=mass, p=power of r
% Entries in the vector of dependent variables are:
% w(1)-position(t), w(2)-velocity(t), w(3)-theta(t)
derivs = [ w(2);-(a/m)*w(1).^p+L^2/m^2./w(1).^3;L/m./w(1).^2];