function [N,H,k,Q,cw,tau]=MODIS_Nd_func_of_re_LWP(W,reff,CTT,CTP)
%function [N,H,k,Q,cw,tau]=MODIS_Nd_func_of_re_LWP(W,reff,CTT,CTP)
%  
% Returns N in cm^{-3} and H in metres as a funciton of W in kg/m2.
% Also returns tau (unitless)
% reff - effective radius from MODIS in metres
% CTT - Cloud Top Temperature (K)
% CTP - Cloud Top Pressure (Pa)
% Created 16th July, 2013 - changed to give as function of W instead of tau
% Daniel Grosvenor (daniel.p.grosvenor@gmail.com)


Q=2;
k=0.8; %Lu&Seinfeld (2006) found it varies between 0.5 and 0.9
        

[cw]=adlwcgm2_just_robs(CTT,CTP) * 1e-3; %*1e-3 to convert to kg/m4

%this value is appropriate for CTT=285K (approx for VOCALS) if want to use a constant value           
%cw=0.002e-3; %kg/m4 - this value is more appropriate - the Fig. 2 in
%Bennartz (2007) is wrong

rhoL=1000; %water density kg/m3

tau = W .*9/5 ./rhoL ./reff;
        %W is LWP (kg/m2)
%        W=5/9*rhoL.*tau.*reff; %this is the one that takes into account adiabaticity
%        W=2/3*rhoL.*tau.*reff; %this one assumes vertical homgeneity


H=sqrt(2*W./cw);  %cloud depth (m)

%droplet number
N = 2.^(-5/2)./k .* tau.^3 .* W.^(-5/2) .* (3*pi.*Q/5).^-3 ...
    .* (3/(4*pi.*rhoL)).^-2 .* sqrt(cw);
%N is proportional to tau.^0.5 * Reff.^(-2.5) - so more sensitive to Reff
%than tau.

N=N/1e6; %convert from m3 to cm3.

