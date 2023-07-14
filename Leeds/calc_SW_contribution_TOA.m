%Calculate SW up TOA for a given cloud Nd, LWP (cloudy sky only) and cloud
%fraction.

%% SW calcluations

function [SW_estimated_forcing_cf_PIPD_TOA_timemean, SW_estimated_forcing_Nd_PIPD_TOA_timemean, SW_estimated_forcing_lwp_PIPD_TOA_timemean, ...
    SW_PI_TOA, SW_PD_TOA] = ...
    calc_SW_contribution_TOA(f0,W0,N0,f1,W1,N1,W0_2,N0_2,W1_2,N1_2,SW_down_PD_ALL,cf_min,A_clear)

i_Liu = 1; %switch to use the Liu droplet size distribution width parameterization for consistency with UKESM.
transmission_atmos = 0.6; %Assume an atmospheric transmissivity - can use this to tune somewhat
transmission_atmos = 0.88; %This value worked best for the ACSIS runs from GC20 for TOA fluxes
transmission_atmos = 0.86; %This value worked best for the ACSIS runs from GC20 for TOA fluxes
transmission_atmos = 0.82; %This value worked best for the ACSIS runs from GC20 for TOA fluxes
transmission_atmos = 0.72; %This value worked best for the ACSIS runs from GC20 for TOA fluxes
transmission_atmos = 0.68; %This value worked best for the ACSIS runs from GC20 for TOA fluxes
transmission_atmos = 0.66; %This value worked best for the ACSIS runs from GC20 for TOA fluxes
%transmission_atmos = 0.95;

SW_in_cloud_top = SW_down_PD_ALL .* transmission_atmos; %downwelling SW at the cloud top (assuming low clouds here)

%PD baseline calc for CF :- SW_estimated_forcing_cf_PD = SW_PD_TOA - SW_cf2;
%[Ac,tau,T_f,A,SW_PI_TOA] = calc_SW(f0,W0,N0,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_PI_TOA] = calc_SW(f0,W0_2,N0_2,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);
%[Ac,tau,T_f,SW_PD_TOA] = calc_SW(f1,W1,N1,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_PD_TOA] = calc_SW(f1,W1_2,N1_2,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);

[Ac,tau,T_f,SW_cf_TOA] = calc_SW(f1,W0_2,N0_2,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_lwp_TOA] = calc_SW(f0,W1,N0,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_Nd_TOA] = calc_SW(f0,W0,N1,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);

[Ac,tau,T_f,SW_cf2_TOA] = calc_SW(f0,W1_2,N1_2,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_lwp2_TOA] = calc_SW(f1,W0,N1,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_Nd2_TOA] = calc_SW(f1,W1,N0,SW_in_cloud_top,A_clear,transmission_atmos,cf_min,0,NaN,i_Liu);



%% Estimated forcings

%Using PI as a baseline
SW_estimated_forcing_cf_TOA = - (SW_cf_TOA - SW_PI_TOA); %minus signs here because an increase in outgoing TOA; -ve for polluted
SW_estimated_forcing_Nd_TOA = - (SW_Nd_TOA - SW_PI_TOA); %is a negative forcing in terms of warming / energy input
SW_estimated_forcing_lwp_TOA = - (SW_lwp_TOA - SW_PI_TOA);

%Using PD as a baseline
SW_estimated_forcing_cf_PD_TOA = - (SW_PD_TOA - SW_cf2_TOA); % - ve for polluted
SW_estimated_forcing_Nd_PD_TOA = - (SW_PD_TOA - SW_Nd2_TOA); 
SW_estimated_forcing_lwp_PD_TOA = - (SW_PD_TOA - SW_lwp2_TOA);

%Take avearges of PI and PD baseline values and time averages
SW_estimated_forcing_cf_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_cf_PD_TOA,3) + meanNoNan(SW_estimated_forcing_cf_TOA,3) )./2;
SW_estimated_forcing_Nd_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_Nd_PD_TOA,3) + meanNoNan(SW_estimated_forcing_Nd_TOA,3) )./2;
SW_estimated_forcing_lwp_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_lwp_PD_TOA,3) + meanNoNan(SW_estimated_forcing_lwp_TOA,3) )./2;


