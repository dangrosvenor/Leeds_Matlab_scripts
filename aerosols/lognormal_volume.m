function dndlnD=lognormal_volume(D,N,sig,Dg)
%function dVlnD=lognormal(D,N,sig,Dg)
%returns dV/dlnD for lognormal distribution of total number N standard
%deviation sig and the median diameter Dg (of the NUMBER distribution)
% E.g. eqn. 8.51 of Seinfeld and Pandis textbook, but multiplied by the
% diameter (D) since 8.51 is dV/dD rather than dV/dlnD
% This follows from eqn. 8.23 which states that dV/dlnD = D * dV/dD
% N.B. Dg is the MEDIAN diameter (rather than mode, mean, etc.)

dndlnD= pi*D.^3.*N./(6*sqrt(2*pi).*log(sig)) * exp( - (log(D/Dg)).^2/2/(log(sig)^2) ); 


