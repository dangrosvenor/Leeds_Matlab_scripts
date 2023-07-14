function N_IN = microphysics_Fletcher_IN(tc)
%function N_IN = IN_Fletcher(tc,P)
%calculates the WRF Fletcher number of IN (function of temperature only)
%tc = temperature in celsius
%N_IN=0.01*exp(-0.6*tc)*1e-3; %per L

N_IN=0.01*exp(-0.6*tc)*1e-3; %per L

