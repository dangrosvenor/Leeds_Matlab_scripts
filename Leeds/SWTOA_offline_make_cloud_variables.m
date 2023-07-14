%Calculate SW up TOA for a given cloud Nd, LWP (cloudy sky only) and cloud
%fraction.

%% Gather and pre-process the cloud variables

function [f0,N0,W0,f1,N1,W1,W0_2,N0_2,W1_2,N1_2] = SWTOA_offline_make_cloud_variables(f0,N0,W0,f1,N1,W1,cf_min)

%it_sw = 1:size(SW_down_PI_ALL,3); %The required time indices

%f0 = low_CF_PI_ALL(:,:,it_sw); f1 = low_CF_PD_ALL(:,:,it_sw);
%N0 = Nd_PI_ALL(:,:,it_sw)/1e6; N1 = Nd_PD_ALL(:,:,it_sw)/1e6; %convert to per cc for calc_SW function
%W0 = LWP_PI_ALL(:,:,it_sw)./f0; W1 = LWP_PD_ALL(:,:,it_sw)./f1; %Need to supply the in-cloud LWP, so divide by low CF
%cf_min = 0.01;
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
W0_2(i) = W1(i); 
N0_2(i) = N1(i); %replace the PI values with the PD ones

%Do a similar thing for when use PD as the baseline - set W and N to PI
%values if going from zero to some CF between PD and PI
i=find(f1<cf_min & f0>=cf_min);
W1_2 = W1; N1_2=N1;
W1_2(i) = W0(i); 
N1_2(i) = N0(i); %replace the PD values with the PI ones



