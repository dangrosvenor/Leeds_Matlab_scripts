function [N1,s]=twomeyconc(w,c,k) 
%N=CS^k relationship from Twomey (1959) - C in cm^-3 and w in m/s

s=(1.63e-3*(w*100)^(1.5)/c/k/beta(1.5,k/2) ).^(1/(k+2));

N1=c*s.^k;

%N2=c.^(2/(k+2)).*(1.63e-3*(w*100)^(3/2)/k/beta(1.5,k/2) ).^(k/(k+2));