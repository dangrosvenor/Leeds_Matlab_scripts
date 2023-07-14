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
Nvals = [25 50 100]; %cm^-3 
Hcb = 500; %Height of cloud base - assuming rho_air*qR is constant up to cloud base - likely not true due to rain evap
 %Fairly large uncertainty here, but can scale linearly.
LWP = [0:10:180]; %g/m2
%Estimate precip rate for various LWPs up to AMSRE threshold for precip

%Set parameters required to calculate the cloud tau
CTT= 280;
rv = 100e-6;
k=0.6; %k for the rain dist
kc = 0.8; %for the cloud distribution
v_drizzle = 2; %m/s  Required as a started guess for the iterative method for all but the Marshall Palmer dist.

cw=MODIS_justN_func(1e9,1e9,'calc',0,CTT,'cw');  %kg/m4 - may be of use.
B = Nd_equations_B_func_of_tau_W(CTT);


for i=1:length(Nvals)
    N=Nvals(i);


    tauC_N50 = ( N .* (LWP/1e3).^(5/2) / B ).^(1/3);
    re_c = re_from_N_LWP(N*1e6,LWP/1e3,kc,CTT); %re in metres


    P = precip_rate_Wood_2008(LWP,N); %returned in mm/hr
    [RWP_01_N50,tau_R_01_N50, re, Vfall, lam] = rwp_from_precip_rate(P,v_drizzle,Hcb,rv,k);
    Wdiff_01_N50 = LWP .* (tau_R_01_N50 ./ tauC_N50);  %Wdiff is the extra LWP from the addition of the rain for MODIS

    LWP_R_modis_01_N50 = 5/9 * 1e3 * re_c .* tau_R_01_N50 * 1e3; %convert to g/m2 -- MODIS LWP from RAINY part
    PB_01_N50 = 100 * (LWP_R_modis_01_N50 - RWP_01_N50) ./ (LWP + RWP_01_N50); %percentage diff LWP MODIS minus AMSRE
    
    PR_save{i} = P;
    RWP_save{i} = RWP_01_N50;
    tauR_save{i} = tau_R_01_N50;
    Wdiff_save{i} = Wdiff_01_N50;
    LWP_R_modis_save{i} = LWP_R_modis_01_N50;
    PB_save{i} = PB_01_N50;
    re_save{i} = re;
    Vfall_save{i} = Vfall;
    lam_save{i} = lam;
    
    labs_RWP(i).l = ['N_d = ' num2str(N) ' cm^{-3}'];

    
end

