function [S,Dcrit,nd,Dspacing,logD]=makeBinsNenes(Dg,sig,N,Nbins,T)
%returns S in e/es form

% Dg=0.02e-6; %mode diameter
% sig=2.5;
% Nbins=200;
% N=200e6; %total no. conc m^-3
% T=298;

Dlow=Dg/10/sig;
Dhigh=Dg*10*sig;
Dspacing=( log(Dhigh)-log(Dlow) )/Nbins; %log spacing

logD=[log(Dlow):Dspacing:log(Dhigh)]; %equally spaced in log space

Dmid=(logD(2:end) + logD(1:end-1) )/2;
Dmid=exp(Dmid);

nd=lognormal(Dmid,N,sig,Dg)*Dspacing; %total no. in each bin

%find critical supersaturation from radius (ammonium sulphate)
[S,Dcrit]=Scrit(T,fliplr(exp(logD)/2)); %fliplr flips matrix so that low supersats will be first
S=S/100; %convert to e/es form for consistency with Nenes.m

logD=fliplr(logD); %flip so that runs in same direction as S