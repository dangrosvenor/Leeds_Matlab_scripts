function [gamma_w,gamma_w2]=moist_ad_lapse_rate(T,P,flag)  
%Moist adiabatic lapse rate in K/m
%from the AMS glossary definition

if nargin==2
    flag='rob';
end

f=1e6*28.97/18; %

%qv=SatVapPress(T,'goff','liq',P,1)/f; %use saturated qv mixing ratio
%qv=qsatw_rob(T,P);

[qv,dq_dT]=SatVapPress2(T,flag,'liq',P,1);
dq_dT = dq_dT/f;
qv=qv/f;

g=9.81;    
R=287.04;    
%Lv=2.5e6;    %Latent heat of vap-orization
Lv=latent_heat_temp_variation(T);
cp=1004.67;    

eps=0.622;


gamma_w = g * (1+Lv.*qv./(R.*T) ) ./ (  cp + (Lv.^2.*qv*eps ./ (R.*T.^2) )  );

%need e and dq_dT


%p-e = eps*P./(qv+eps)

%gamma_w2 = g * (1+Lv.*qv.*P/(R.*T.*eps.*P./(qv+eps) ) ) ./ (  cp + (Lv.*dq_dT) );  =
gamma_w2 = g * ( 1 + Lv.*qv.*(qv+eps)./(R.*T.*eps) ) ./ (  cp + (Lv.*dq_dT) );


