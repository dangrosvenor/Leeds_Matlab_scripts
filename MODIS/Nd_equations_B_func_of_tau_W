function re = re_for_precip_rate(tau,P,CTT)
% function re = re_for_precip_rate(tau,P,CTT)
% Returns re in microns for a given precip rate and tau using the formula
% for precip rate given in Wood (2008).
% P is the precip rate in mm/hr and CTT (optional) is the cloud top
% temperature
% Are assuming that LWP = 5/9*rhoW*tau*re here (i.e. the adiabatic
% assumption rather then vertically homogeneous).

if ~exist('CTT')
    CTT=280; %pick a sensible value
end


%Calculate the factor for the Nd as a function of tau and re realtionship
%This a function of CTT only for this calc (P is assumed in MODIS_justN_func
re_set = 14e-6; %can pick any and tau value.
N=MODIS_justN_func(tau,re_set,'calc',0,CTT,'N');
K=N/(tau.^0.5*re_set.^-2.5);


re = 1e6*( (P*24).^0.1633 / tau.^0.1429 *  ( K / (1e6*0.37*5/9) ).^0.2857 );