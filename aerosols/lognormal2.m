function dndlnD=lognormal(logD,N,sig,Dg)
%returns dn/dlnD for lognormal distribution of total number N standard deviation sig and number mode diameter Dg
D=exp(logD);
dndlnD= N/sqrt(2*pi)/log(sig) * exp( - (log(D/Dg)).^2/2/(log(sig)^2) );