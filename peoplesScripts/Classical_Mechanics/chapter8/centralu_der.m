%centralu_der.m: returns the derivatives for 
%-a*r^p central force, where r=1/u
function derivs = centralu_der( h, v, flag,a,L,m,p)
% a=force strength, L=angular momentum, m=mass, p=power of r
% Entries in the vector of dependent variables are:
% v(1) is u(h), v(2) is u'(t), here h=angle, u is a function of angle
derivs = [v(2);-v(1)-((m/L^2)./v(1)^2).*(-a./v(1)^p)];