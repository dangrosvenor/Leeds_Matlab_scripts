function [dqldz]=adlwcgm2_just_robs(t,p)
%gives condensation rate dqldz in g/m3 per m
%[dqldz]=adlwcgm2_just_robs(t,p)
%p in Pa
%t in K
  
qs=qsatw_rob(t,p); %water saturation vapour pressure kg/kg

f=1.6094e+06;       
rho=p./(287.04*t);    
g=9.81;    
R=287.04;    
L=2.5e6;    
cp=1004.67;    

eps=0.622;
      
dqldz=(g*qs./(R.*t)).*(L*eps./(cp*t)-1.0)./(1+eps.*L.*L.*qs./(R.*t.*t.*cp));   
dqldz=rho.*dqldz.*1000.; %dqldz in g/m3 per m

  

