function dndlnD=lognormal(D,N,sig,Dg)
%function dndlnD=lognormal(D,N,sig,Dg)
%returns dn/dlnD for lognormal distribution of total number N standard deviation sig and the median diameter Dg
% E.g. eqn. 8.33 of Seinfeld and Pandis textbook.
% N.B. Dg is the MEDIAN diameter (rather than mode, mean, etc.)

dndlnD= N/sqrt(2*pi)/log(sig) * exp( - (log(D/Dg)).^2/2/(log(sig)^2) );

%Also, N.B. dn/dD is given by :-

% dndD= N./sqrt(2*pi)./log(sig)./D .* exp( - (log(D/Dg)).^2/2/(log(sig)^2) );

% I.e. the same, but divided by D