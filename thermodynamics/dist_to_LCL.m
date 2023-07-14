function [dist_to_CB,T_CB] = dist_to_LCL(T,P,RH)

f = 1.6094e+06;
qsat = SatVapPress(T,'goff','liq',P,1)/f; %kg/kg
q_parcel = RH/100 * qsat;
T_CB=Tdew(q_parcel,P); %the dew point temperature for our q
LR_dry = 9.81/1005;  % K/m
dist_to_CB = (T - T_CB)/LR_dry; %in metres
