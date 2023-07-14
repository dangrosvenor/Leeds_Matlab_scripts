function th = potemp(T,P)
%function th = potemp(T,P)

th = T.*( (1000e2./P).^0.286 );