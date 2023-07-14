function [W]=MODIS_justN_func(tau,reff,Wflag,WMOD,CTT,out_flag)
%function [W]=MODIS_justN_func(tau,reff,Wflag,WMOD,CTT,out_flag)
%set what to output in out_flag - 'N', 'H' or 'W'
%Returns N in cm^{-3} and H in metres
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

Q=2;
k=0.8; %Lu&Seinfeld (2006) found it varies between 0.5 and 0.9
        
P=850*1e2; %Convert to Pa. Using constant pressure here. Could vary according to the pressure also? Although
%the pressure dependence is not very strong
Parr=ones(size(CTT))*P; %

switch out_flag
    case 'cw'
        %just output this and exit
        [W]=adlwcgm2_just_robs(CTT,Parr) * 1e-3; %*1e-3 to convert to kg/m4     
        return
    otherwise
        [cw]=adlwcgm2_just_robs(CTT,Parr) * 1e-3; %*1e-3 to convert to kg/m4
end

%this value is appropriate for CTT=285K (approx for VOCALS) if want to use a constant value           
%cw=0.002e-3; %kg/m4 - this value is more appropriate - the Fig. 2 in
%Bennartz (2007) is wrong

rhoL=1000; %water density kg/m3

switch Wflag
    case 'calc'

        %W is LWP (kg/m2)
        W=5/9*rhoL.*tau.*reff; %this is the one that takes into account adiabaticity
%        W=2/3*rhoL.*tau.*reff; %this one assumes vertical homgeneity

        %Rob says that don't need to scale W by cf if are using the MODIS W as this
        %is already scaled - i.e. is the in-cloud rather than grid-box
        %average

    case 'MODIS'
        W=WMOD; %use the MODIS supplied value (already scaled for cf)
end

switch out_flag
    case 'H'

        W=sqrt(2*W./cw);  %cloud depth (m)

    case 'N'

        %droplet number (cm3)
        W = 1e-6 * 2.^(-5/2)./k .* tau.^3 .* W.^(-5/2) .* (3*pi.*Q/5).^-3 ...
            .* (3/(4*pi.*rhoL)).^-2 .* sqrt(cw);
        %N is proportional to tau.^0.5 * Reff.^(-2.5) - so more sensitive to Reff
        %than tau.


end

