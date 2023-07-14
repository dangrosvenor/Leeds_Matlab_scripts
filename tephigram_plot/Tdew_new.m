function y=Tdew_zero_function(T,Q,P)        
%function y=Tdew_zero_function(T,Q,P)    
%Q is the vapour MR in kg/kg
%P is the air pressure in Pa
%T is the temperature in K

f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
y=SatVapPress(T,'goff','liq',P,1)/f - Q;

