function N_IN = microphysics_WRF_IN(tc)
%function N_IN = IN_wrf(tc,P)
%calculates the WRF Cooper number of IN (function of temperature only)
%tc = temperature in celsius
%N_IN = 0.005*exp(0.304*(-tc)); %per L

N_IN = 0.005*exp(0.304*(-tc)); %per L

