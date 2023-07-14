%Calculate SW up TOA for a given cloud Nd, LWP (cloudy sky only) and cloud
%fraction.

%% Gather and pre-process the cloud variables

function calc_cloud_variables_for_SW_TOA_calc()

it_sw = 1:size(SW_down_PI_ALL,3); %The required time indices

f0 = low_CF_PI_ALL(:,:,it_sw); f1 = low_CF_PD_ALL(:,:,it_sw);
N0 = Nd_PI_ALL(:,:,it_sw)/1e6; N1 = Nd_PD_ALL(:,:,it_sw)/1e6; %convert to per cc for calc_SW function
W0 = LWP_PI_ALL(:,:,it_sw)./f0; W1 = LWP_PD_ALL(:,:,it_sw)./f1; %Need to supply the in-cloud LWP, so divide by low CF
cf_min = 0.01;
W0(f0<cf_min)=NaN; W1(f1<cf_min)=NaN; %regions with no cloud - set W to NaN since have no cf to divide by
N0(f0<cf_min)=NaN; N1(f1<cf_min)=NaN; %Do the same for Nd for consistency

%But only NaN the cloud fraction when we are not at zero CF - for zero CF W and N are allowed
%to be NaN.
i=find(isnan(N0)==1 & f0>=cf_min); W0(i)=NaN; 
i=find(isnan(W0)==1 & f0>=cf_min); N0(i)=NaN; 
i=find(isnan(N1)==1 & f1>=cf_min); W1(i)=NaN; 
i=find(isnan(W1)==1 & f1>=cf_min); N1(i)=NaN; 

%Special case where we have no CF in PI, but some cloud in PD - this would
%usually lead to the change in CF being ignored since the PD W and N would
%be NaN. So, set the PI W and N values in these cases to be equal to the
%PD values, so that this effect is incorporated into the CF change effect
i=find(f0<cf_min & f1>=cf_min); %If they are both <cf_min then they will both be treated as clear-sky in calc_SW
W0_2 = W0; N0_2=N0;
W0_2(i) = W1(i); N0_2(i) = N1(i); %replace the PI values with the PD ones

%Do a similar thing for when use PD as the baseline - set W and N to PI
%values if going from zero to some CF between PD and PI
i=find(f1<cf_min & f0>=cf_min);
W1_2 = W1; N1_2=N1;
W1_2(i) = W0(i); N1_2(i) = N0(i); %replace the PD values with the PI ones





%% SW calcluations

[SW_estimated_forcing_cf_PIPD_TOA_timemean, SW_estimated_forcing_Nd_PIPD_TOA_timemean, SW_estimated_forcing_lwp_PIPD_TOA_timemean] = ...
    function calc_SW_contribution_TOA(f0,W0,N0,SW_down_PD_ALL,cf_min)

i_Liu = 1; %switch to use the Liu droplet size distribution width parameterization for consistency with UKESM.
transmission_atmos = 0.6; %Assume an atmospheric transmissivity - can use this to tune somewhat
transmission_atmos = 0.88; %This value worked best for the ACSIS runs from GC20

SW_in_TOA = SW_down_PD_ALL .* transmission_atmos; %downwelling SW at the cloud top (assuming low clouds here)

%PD baseline calc for CF :- SW_estimated_forcing_cf_PD = SW_PD_TOA - SW_cf2;
[Ac,tau,T_f,SW_PI_TOA] = calc_SW(f0,W0,N0,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_PD_TOA] = calc_SW(f1,W1,N1,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);

[Ac,tau,T_f,SW_cf_TOA] = calc_SW(f1,W0_2,N0_2,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_lwp_TOA] = calc_SW(f0,W1,N0,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_Nd_TOA] = calc_SW(f0,W0,N1,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);

[Ac,tau,T_f,SW_cf2_TOA] = calc_SW(f0,W1_2,N1_2,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_lwp2_TOA] = calc_SW(f1,W0,N1,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);
[Ac,tau,T_f,SW_Nd2_TOA] = calc_SW(f1,W1,N0,SW_in_TOA,transmission_atmos,cf_min,0,NaN,i_Liu);



%% Estimated forcings

%Using PI as a baseline
SW_estimated_forcing_cf_TOA = SW_cf_TOA - SW_PI_TOA;
SW_estimated_forcing_Nd_TOA = SW_Nd_TOA - SW_PI_TOA;
SW_estimated_forcing_lwp_TOA = SW_lwp_TOA - SW_PI_TOA;

%Using PD as a baseline
SW_estimated_forcing_cf_PD_TOA = SW_PD_TOA - SW_cf2_TOA;
SW_estimated_forcing_Nd_PD_TOA = SW_PD_TOA - SW_Nd2_TOA; 
SW_estimated_forcing_lwp_PD_TOA = SW_PD_TOA - SW_lwp2_TOA;

%Take avearges of PI and PD baseline values and time averages
SW_estimated_forcing_cf_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_cf_TOA_PD,3) + meanNoNan(SW_estimated_forcing_cf_TOA,3) )./2;
SW_estimated_forcing_Nd_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_Nd_TOA_PD,3) + meanNoNan(SW_estimated_forcing_Nd_TOA,3) )./2;
SW_estimated_forcing_lwp_PIPD_TOA_timemean = ( meanNoNan(SW_estimated_forcing_lwp_TOA_PD,3) + meanNoNan(SW_estimated_forcing_lwp_TOA,3) )./2;


