function [N,H,W,k,Q,cw]=MODIS_N_H_func(tau,reff,Wflag,WMOD,CTT,fad,k)
%[N,H,W,k,Q,cw]=MODIS_N_H_func(tau,reff,Wflag,WMOD,CTT)
%Returns N in cm^{-3} and H in metres, W in kg/m2
%cf = cloud fraction - Have realised that this isn't needed after speaking
%to Rob. Bennartz (2007) used it in his equation, but this was only for the
%AMSR microwave satellite. MODIS is already scaled for cloud fraction.
%CTT = cloud top temperature for calculation of cw
%tau - optical depth from MODIS (no units)
%reff - effective radius from MODIS in metres
%Wflag - flag to say whether to calculate the LWP (=W) from tau and reff (Wflag='calc' or
%to use the MODIS value (Wflag='MODIS')
%
%updated 27th Oct 2011 - now includes the dependence of the condensation
%rate (cw) on temperature - reduces significantly with lowering temperatures

if ~exist('fad')
    fad=1.0;
end

Q=2;

if ~exist('k') | isnan(k)==1
    k=0.8; %Lu&Seinfeld (2006) found it varies between 0.5 and 0.9
end

P=850*1e2; %Convert to Pa. Using constant pressure here. Could vary according to the pressure also? Although
%the pressure dependence is not very strong
Parr=ones(size(CTT))*P; %
[cw]=adlwcgm2_just_robs(CTT,Parr) * 1e-3; %*1e-3 to convert to kg/m4

%this value is appropriate for CTT=285K (approx for VOCALS) if want to use a constant value           
%cw=0.002e-3; %kg/m4 - this value is more appropriate - the Fig. 2 in
%Bennartz (2007) is wrong

rhoL=1000; %water density kg/m3

switch Wflag    
    case 'MODIS'
        W=WMOD; %use the MODIS supplied value (already scaled for cf)
        
    otherwise

        %W is LWP (kg/m2)
        W=5/9*rhoL.*tau.*reff; %this is the one that takes into account adiabaticity
%        W=2/3*rhoL.*tau.*reff; %this one assumes vertical homgeneity

        %Rob says that don't need to scale W by cf if are using the MODIS W as this
        %is already scaled - i.e. is the in-cloud average, not the gridbox
        %average     
end

H=sqrt(2*W./cw);  %cloud depth (m)

%droplet number
N = 2.^(-5/2)./k .* tau.^3 .* W.^(-5/2) .* (3*pi.*Q/5).^-3 ...
    .* (3/(4*pi.*rhoL)).^-2 .* sqrt(fad*cw);
%N is proportional to tau.^0.5 * Reff.^(-2.5) - so more sensitive to Reff
%than tau.

%N.B. - above gives the same answer to that in the final version of the
%2018 review paper (Grosvenor Rev Geo Phys) when use tau and reff :-
%N = sqrt(5)./(2.*pi.*k) .* (fad.*cw.*tau ./ (Q.*rhoL.*reff.^5)).^0.5;

N=N/1e6; %convert from m3 to cm3.

