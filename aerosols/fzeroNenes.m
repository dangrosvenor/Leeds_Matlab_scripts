function I=fzeroNenes(smax,a)
%use to run with fzero e.g. fzero(@fzeronenes,Smax_guess)
%see multiSmaxs for routine for plotting viewing different distributions

Dg=0.02e-6; %mode diameter (m)
sig=2.5;
Nbins=200;
N=200e6; %total no. conc m^-3
T=298;

Dlow=Dg/10/sig;
Dhigh=Dg*10*sig;
Dspacing=( log(Dhigh)-log(Dlow) )/Nbins; %log spacing

logD=[log(Dlow):Dspacing:log(Dhigh)]; %equally spaced in log space

Dmid=(logD(2:end) + logD(1:end-1) )/2;
Dmid=exp(Dmid);

nd=lognormal(Dmid,N,sig,Dg); %total no. in each bin


[S,Dcrit]=Scrit(T,fliplr(exp(logD)/2)); %fliplr flips matrix so that low supersats will be first


[I,sp,smax]=nenes(smax,10e4,283,10,nd*Dspacing,S);
