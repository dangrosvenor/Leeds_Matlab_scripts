function [dqldz,dqldz_dan,dLWCdT,dqldz_dan2,dqldz_dan3,dqldz_dan4,dqdz_num,dqldz_dan5,dqldz_dan6,dqldz_dan7,dqldz_dan_goff,dqldz_Lc,flag,dqldz_sami]=adlwcgm2_rob(t,p)
%p in Pa
%t in K
f=1.6094e+06;

%qs=qsatw_rob(t,p); 
qs=SatVapPress(t,'Rob','liq',p,1)/f;
[e,de_dT]=SatVapPress2(t,'Rob','liq');
[moist_ad_1,moist_ad_2] = moist_ad_lapse_rate(t,p);
%moist_ad_2 is the moist adiabatic temperature without the
%C-C approximation. Also using temperature varation of L.
%moist_ad_1 uses the C-C relationship (with T variation of L)


%Clausius-Clapeyron
%Is actually pretty close to the numerical value from differentiation of qs
% IF the temperature variation of L is taken into account.
%de_dT_CC=L*eps .*e /R ./t.^2;

      
rho=p./(287.04*t);    
g=9.81;    
R=287.04;    
%L=2.5e6;  
L=latent_heat_temp_variation(t);
cp=1004.67;    

eps=0.622;

%gam in g/kg per metre    
    
%gam in g/m3/m    
    

dqldz=(g*qs./(R.*t)).*(L*eps./(cp*t)-1.0)./(1+eps.*L.*L.*qs./(R.*t.*t.*cp));
dqldz=rho.*dqldz.*1000.;

%from Sami's manuscript (this is wrong)
dqldz_sami = 1e3 * 1./(R*t/eps) .* de_dT .* moist_ad_2;  %g/m3/m
%comes from the assumption that dq/dz = 1./(R*t/eps) * de/dz  
%(from q(kg/kg) = eps*e/P, then he assumed dq/dz = eps/P * de/dz 
% * rho (=p/RT) to turn into kg/m3 gives dql/dz(kg/m4) = eps/RT * de/dz as
% above - BUT since q is a function of P as well as T (from e=e(T)) you
% can't do this). Need to use the total derivative (see Condensation rate Albrecth document for this 
%type of derivation).

%dqldz in g/m3 per m  

Lc=2.5e6;

dqldz_Lc=(g*qs./(R.*t)).*(Lc*eps./(cp*t)-1.0)./(1+eps.*Lc.*Lc.*qs./(R.*t.*t.*cp));
dqldz_Lc=rho.*dqldz_Lc.*1000.;


%numerical integration of qs(T,P)
%for a step dz calculate dT from moist ad lapse rate dT/dz
%also calculate dp from dp/dz=-rho*g (hydrostatic equation)
%then use the new T and P to calculate the change in qs
%which can be divided by dz to give dq/dz

dz=1; %small increase in z in metres

for i=1:length(t)
    dP = -rho(i)*g*dz;
    dT = -moist_ad_lapse_rate(t(i),p(i))*dz;
    qnew = qsatw_rob(t(i)+dT,p(i)+dP); %use new p and T
    qold = qsatw_rob(t(i),p(i)); %use old p and T
%    qnew = SatVapPress(t(i)+dT,'goff','liq',p(i)+dP,1)/f; %use new p and T    
%    qold = SatVapPress(t(i),'goff','liq',p(i),1)/f; %use new p and T    
    dq = qold - qnew; 
    dqdz_num(i) = dq/dz;
    
end

dqdz_num = dqdz_num.*rho*1000; %convert to g/m4




dp=0.5; %hPa

for i=1:length(p)
    %test using another approx method
    %first estimate the gradient of adiabatic LWC with temperature
    %this method matches the correct one (above), but only at colder
    %temperatures (around 290K)??
    
    p2=[p(i) p(i)-dp]; %calculate at a cloudy point and at a location just above
    %could be at any point in the cloud as we just want the gradient of LWC vs. T
    %TT
    [LWC,Tad]=adLWC_PaulLawson_simple(p2(1),t(i),p2);
    dLWCdT(i) = -( LWC(2) - LWC(1) ) / ( Tad(2) - Tad(1) ); %gradient - minus sign since LWC goes up
    %as temperautre goes down, in g/m3/K

    %so dLWC/dz should then be dLWCdT * dT/dz
    %approx moist lapse rate from AMS
    gw=moist_ad_lapse_rate(t(i),p(i)); %K/m

    dqldz_dan(i) = dLWCdT(i)*gw;

end

%another attempt
%from p.8 of Byers assuming that dw/dz=dw/dp * dp/dz, substituting in
%hydrostatic equation and assuming dT/dz = dT/dp*dp/dz
%this method works and can be shown to be the same as Rob's method (sub in
%moist ad equation)

dqldz_dan2 = 1./L.* (g - cp*moist_ad_1 );
dqldz_dan2 = dqldz_dan2 .* rho *1000; %convert to g/m4 from kg/kg/m

%yet another attempt! - this didn't seem to work
%from p.31 of Byers assuming that using Clasius-Claperyron equation
%de/dT=Le/RT^2. Converted to q using q=e*eps/(P-e) leads to dq/dT=...
dqdT = L.*qs.*(qs+eps) ./ (R.*t.^2.*eps);  %(kg/kg/K)
%then use dq/dz = dq/dT * dT/dz

dqldz_dan3 = dqdT.*moist_ad_1*1000.*rho;
%convert to kg/m4 by the 1000*rho factor



%2nd Feb 2012 - trying this again with new derivation of total derivative
%using total derivative since q=q(T,P), a funtion of pressure and
%temperature. this is the same as Rob's at lower temperatures, but then becomes
%different. I think it is because
%the clausius_clapeyron relation is inexact. It approximates the
%gradient of e(T) - but can be solved analytically and so negates the need
%for the complicated qsat formlulae. See thermodynamics/clausius_clapeyron.m
%gets the same result as in Albrecht paper

%Albrecht's formula
%this is dqv/dz
%now using moist_ad_2 - i.e. the moist adiabatic temperature without the
%C-C approximation. Also using temperature varation of L
dqldz_dan5 = qs.*p.*g./(R.*t.*(p-e)) - qs.*(qs+eps).*L.*moist_ad_2./(R.*t.^2) ;

%dq_L/dz is negative dqv/dz
dqldz_dan5 = - dqldz_dan5 * 1000 .*rho;

% Could use the total derivative method above, but calculate
% de/dT instead of using C-C
% dq_dT = ( (qs+eps)./(p-e) ) * de/dT
% dq_dp = -qs ./ (p-e)
% then dq/dz = dq/dp * dp/dz + dq/dT*dT/dz



% dT=0.01;
% e = qs.*p ./ (qs + eps);
% qs_dT=qsatw_rob(t+dT,p); 
% e_dT = qs_dT.*p ./ (qs_dT + eps);
% de_dT = (e_dT - e) / dT;


%Albrecht's total derivative method, but using numerical de/dT instead of C-C
dqldz_dan6 = qs./(p-e) .*p.*g ./(R*t) - ( (qs+eps)./(p-e) ) .* de_dT .* moist_ad_1;
dqldz_dan6 = - dqldz_dan6 * 1000 .*rho;
%Note, the derivation of the moist adiabat also uses the C-C
%approximation - should be possible to make a version that uses the q
%formula to get dq/dT

%same, but with numerical (or rather based on qs formulae) moist adiabat
%moist_ad_2 is the moist adiabat without the C-C approximation
dqldz_dan7 = qs./(p-e) .*p.*g ./(R*t) - ( (qs+eps)./(p-e) ) .* de_dT .* moist_ad_2;
dqldz_dan7 = - dqldz_dan7 * 1000 .*rho;


dqldz_dan4=0; %keep this to keep the outputs the same

%now test different qs formulations
flag='goff';
%flag='buck2';
%flag='wmo';
%flag='roger';
qs=SatVapPress(t,flag,'liq',p,1)/f;
[e,de_dT]=SatVapPress2(t,flag,'liq');
[moist_ad_1,moist_ad_2] = moist_ad_lapse_rate(t,p,flag);

%moist_ad_2 is the moist adiabat without the C-C approximation
dqldz_dan_goff = qs./(p-e) .*p.*g ./(R*t) - ( (qs+eps)./(p-e) ) .* de_dT .* moist_ad_2;
dqldz_dan_goff = - dqldz_dan_goff * 1000 .*rho;












      
    

