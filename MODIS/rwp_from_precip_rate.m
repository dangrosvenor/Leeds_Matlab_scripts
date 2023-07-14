function [RWP, tauR, re, V2, lam] = rwp_from_precip_rate(P,v_drizzle,Hcb,rv,k)
%Estimates of RWP for given precip rates
%Based on formula for precip flux of P (kg/m2/s) = v_drizzle * rho_air * qR
%Note that kg/m2 is the same as mm w.e. since 1 mm of water over 1 m2
%weighs 1e-3*rho_water = 1 kg
%where v_drizzle is the fall speed of the rain, qR is the rain mixing ratio
%(kg/kg).
%Comstock (2004) shows graph from Pruppacher and Klett of fall speed vs
%drizzle drop radius - v likely to be 0.2-0.3 m/s.
%This allows calc of rho_air*qR. If assume this is constant with height
%(may not be due to rain evap) then can integrate to get RWP. Need to
%assume a cloud base height too.
%
%P is in mm/hr, v_drizzle in m/s and Hcb in metres. 
%RWP returned in g/m2
%Also returns the optical depth (tau) of the rainy part
%This requires an input of rv, the mean droplet radius of the rain in
%metres. And also k = (rv/re)^3.

method='gamma'; % (inc. Marshall Palmer)
%method='Comstock';

if ~exist('rv')
    rv = 60e-6;
end

if ~exist('k')
    k = 0.6;
end

rho_air = 1;
rho_water = 1e3;


switch method
    case 'old'
        qR_vol=P/3600 / (v_drizzle); %converts P to kg/m2/s and then divides by v to get qR*rho_air (kg/m3)
    case 'gamma'
        [V2,qR_vol,nit,re,lam]=fall_speed_V_qR_iteration(P,v_drizzle);  %mass weighted fall speed, iterative calc
    case 'Comstock'
        rmean=rv;  %This is assumed for now
        r0=20e-6;
        M3=Moments_for_wood2005_dist(3,rmean,20e-6);
        Md=Moments_for_wood2005_dist(1.4+3,rmean,20e-6);        
        
        qR_vol = P/3600 / 2.2e5 .* M3./Md;
        
        M2=Moments_for_wood2005_dist(2,rmean,20e-6);
        
        re = M3 ./ M2; %re is the ratio of 3rd to second moments
        
        V2 = P/3600 ./qR_vol;
        lam=0; %lambda not defined for this type of distribution.
        
end

RWP = qR_vol * Hcb * 1e3;  %integrate this with height and multiply by 1e3 to get RWP in g/m2

%tau  = 3/2  / rho_water * k.^(1/3) .*qR_vol ./rv .* Hcb;  %have cancelled out rho_air here since qR_vol is rho_air*qR

%Now work out re from the assumed raindrop distribution
tauR  = 3/2  / rho_water .*qR_vol ./re .* Hcb;  %have cancelled out rho_air here since qR_vol is rho_air*qR
