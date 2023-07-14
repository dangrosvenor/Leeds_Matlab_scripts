function [qsw,esw]=qsatw(t,p)

%output kg/kg
%Saturated vap press wrt water.
%use fns given in Unified Model Documentation No 29
%Calculation of saturated specific humidity and large scale cloud.
%RNB Smith et al. 1990
%t is temperature in K
%p is pressure in Pa
%esw is the saturated vapour pressure in Pa
%qsw is the saturated specific humidity wrt water in kg/kg

t0=273.16;

log10esw=10.79574*(1-t0./t)-5.028*log10(t./t0)+1.50475e-4.*(1-10.^(-8.2369.*(t./t0-1)))+0.42873e-3.*(10.^(4.76955.*(1-t0./t))-1)+2.78614;

esw=10.^(log10esw);
qsw=0.62198.*esw./(p-esw);

