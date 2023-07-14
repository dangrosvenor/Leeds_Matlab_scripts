%Estimates of RWP for given precip rates
%Based on formula for precip flux of P (kg/m2/s) = v_drizzle * rho_air * qR
%where v_drizzle is the mass-weighted fall speed of the rain, qR is the rain mixing ratio
%(kg/kg).
%Can calculate v_drizzle from the assumed dist, fall speed relation and P
%This allows calc of rho_air*qR. If assume this is constant with height
%(may not be due to rain evap) then can integrate to get RWP. Need to
%assume a cloud base height too.

iplot_rwc=0;

%Will plot vs LWP for different fall speed value and N values
N = 50; %cm^-3  - set Nd to 50 per cc for convienience
Hcb = 500; %Height of cloud base - assuming rho_air*qR is constant up to cloud base - likely not true due to rain evap
 %Fairly large uncertainty here, but can scale linearly.
LWP = [0:10:180]; %g/m2
%Estimate precip rate for various LWPs up to AMSRE threshold for precip

%Set parameters required to calculate the cloud tau
CTT= 280;
rv = 100e-6;
k=0.6; %k for the rain dist
kc = 0.8; %for the cloud distribution

cw=MODIS_justN_func(1e9,1e9,'calc',0,CTT,'cw');  %kg/m4

B = Nd_equations_B_func_of_tau_W(CTT);
tauC_N50 = ( N .* (LWP/1e3).^(5/2) / B ).^(1/3);

re_c = re_from_N_LWP(N*1e6,LWP/1e3,kc,CTT); %re in metres

v_drizzle = 2; %m/s
P = precip_rate_Wood_2008(LWP,N); %returned in mm/hr
[RWP_01_N50,tau_R_01_N50] = rwp_from_precip_rate(P,v_drizzle,Hcb,rv,k);
Wdiff_01_N50 = LWP .* (tau_R_01_N50 ./ tauC_N50);  %Wdiff is the extra LWP from the addition of the rain for MODIS

LWP_R_modis_01_N50 = 5/9 * 1e3 * re_c .* tau_R_01_N50 * 1e3; %convert to g/m2 -- MODIS LWP from RAINY part
PB_01_N50 = 100 * (LWP_R_modis_01_N50 - RWP_01_N50) ./ (LWP + RWP_01_N50); %percentage diff LWP MODIS minus AMSRE


v_drizzle = 3; %m/s
P = precip_rate_Wood_2008(LWP,N);
[RWP_02_N50,tau_R_02_N50] = rwp_from_precip_rate(P,v_drizzle,Hcb,rv,k);
Wdiff_02_N50 = LWP .* (tau_R_02_N50 ./ tauC_N50);

LWP_R_modis_02_N50 = 5/9 * 1e3 * re_c .* tau_R_02_N50 * 1e3; %convert to g/m2 -- MODIS LWP from RAINY part
PB_02_N50 = 100 * (LWP_R_modis_02_N50 - RWP_02_N50) ./ (LWP + RWP_02_N50);

if iplot_rwc==1

figure
plot(LWP,RWP_01_N50,'k-');
hold on
plot(LWP,RWP_02_N50,'k--');



xlabel('LWP (g m^{-2})');
ylabel('RWP (g m^{-2})');

end

Nd_str01 = ['N_d = ' num2str(N) ' cm^{-3}'];


%Will plot vs LWP for different fall speed value and N values
N = 25; %cm^-3  - set Nd to 50 per cc for convienience
Hcb = 500; %Height of cloud base - assuming rho_air*qR is constant up to cloud base - likely not true due to rain evap
 %Fairly large uncertainty here, but can scale linearly.


v_drizzle = 2; %m/s
P = precip_rate_Wood_2008(LWP,N);


re_c = re_from_N_LWP(N*1e6,LWP/1e3,kc,CTT); %re in metres

%RWP and tau for the rainy part:-
[RWP_01_N100, tau_R_01_N100] = rwp_from_precip_rate(P,v_drizzle,Hcb,rv,k);
%Tau for the cloudy part only:-
tauC_N100 = ( N .* (LWP/1e3).^(5/2) / B ).^(1/3);  %this is the same regardless of v_drizzle
Wdiff_01_N100 = LWP .* (tau_R_01_N100 ./ tauC_N100);

LWP_R_modis_01_N100 = 5/9 * 1e3 * re_c .* tau_R_01_N100 * 1e3;
PB_01_N100 = 100 * (LWP_R_modis_01_N100 - RWP_01_N100) ./ (LWP + RWP_01_N100);


v_drizzle = 3; %m/s
P = precip_rate_Wood_2008(LWP,N);
[RWP_02_N100, tau_R_02_N100, re, Vfall, lam] = rwp_from_precip_rate(P,v_drizzle,Hcb,rv,k);
Wdiff_02_N100 = LWP .* (tau_R_02_N100 ./ tauC_N100);

LWP_R_modis_02_N100 = 5/9 * 1e3 * re_c .* tau_R_02_N100 * 1e3;
PB_02_N100 = 100 * (LWP_R_modis_02_N100 - RWP_02_N100) ./ (LWP + RWP_02_N100);

Nd_str02 = ['N_d = ' num2str(N) ' cm^{-3}'];
%labs_RWP(1).l=['v=0.2 m/s ' Nd_str01];
labs_RWP(1).l=[Nd_str01];
%labs_RWP(2).l=['v=0.3 m/s ' Nd_str01];
%labs_RWP(2).l=['v=0.2 m/s ' Nd_str02];
labs_RWP(2).l=[Nd_str02];
%labs_RWP(4).l=['v=0.3 m/s ' Nd_str02];

leg_cell = {labs_RWP(1).l, labs_RWP(2).l, labs_RWP(3).l, labs_RWP(4).l};

if iplot_rwc==1
    %figure
    plot(LWP,RWP_01_N100,'b-');
    %hold on
    plot(LWP,RWP_02_N100,'b--');


    %xlabel('LWP (g m^{-2})');
    %ylabel('RWP (g m^{-2})');



    legend(leg_cell,'location','NorthWest');




    figure
    plot(LWP,-100*(RWP_01_N50 - Wdiff_01_N50)./LWP, 'k-');
    hold on
    plot(LWP,-100*(RWP_02_N50 - Wdiff_02_N50) ./LWP, 'k--');
    plot(LWP,-100*(RWP_01_N100 - Wdiff_01_N100) ./LWP, 'b-');
    plot(LWP,-100*(RWP_02_N100 - Wdiff_02_N100) ./LWP, 'b--');



    xlabel('LWP (g m^{-2})');
    %ylabel('LWP difference AMSRE minus MODIS (g m^{-2})');
    ylabel('MODIS bias (%)');

    legend(leg_cell,'location','SouthWest');


end

