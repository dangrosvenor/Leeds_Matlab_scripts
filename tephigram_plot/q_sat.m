function qs=q_sat(T,P)
%saturation mixing ratio in kg/kg
%T is the temperature in K, P the air pressure in Pa
qs=0.622.*GGEW(T)./(P-GGEW(T));
%GGEW calculates the vapour pressure