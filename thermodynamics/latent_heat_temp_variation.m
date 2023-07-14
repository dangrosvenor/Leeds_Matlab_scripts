function LL=latent_heat_temp_variation(T)
%T in K, LL in J/kg
%from Wikipedia - apparently a cubic fit from Rogers and Yau, Short Course
%on cloud physics. For -40 to 40 degC temperature range. Varies from approx
%2600-2400 over this temp range


T=T-273.15;
LL=-0.0000614342.*T.^3 + 0.00158927.*T.^2 - 2.36418.*T + 2500.79;
LL=LL*1e3;