function T=th2t(p,th)
%function T=th2t(p,th)
%calculates temperature from the potemp and pressure:- T=th./( (1000e2./p).^0.286 );

T=th./( (1000e2./p).^0.286 );