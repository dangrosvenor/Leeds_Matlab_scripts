function [p_i,gradient]=SatVapPress2(T,flag,flag2,P,ppmv);
%gives saturation vapour pressure in Pa
%usage: SatVapPress(T,flag,flag2,P,ppmv); P in Pa, T in K
%flag= e.g. 'goff', flag2= 'ice' or 'liq'
%pressure only required if want sat mixing ratio in ppmv
%then set ppmv=1 (optional flag)
%available methods:
%ice: 'buck2' 'buck' 'goff' 'marti' 'teten' 'hyland' 'murphy' case 'lem'
%liq: 'goff' 'bolton' 'roger''buck2' 'buck1' 'wmo' 'hyland' 'sonntag'
%'teten', 'mesonh
%and also the gradient dq/dT

dT = 0.01;

if nargin==3
  
p_i = SatVapPress(T,flag,flag2);

%now calculate the gradient

val_a = SatVapPress(T-dT,flag,flag2);
val_b = SatVapPress(T+dT,flag,flag2);  
    
else
    

p_i = SatVapPress(T,flag,flag2,P,ppmv);

%now calculate the gradient

val_a = SatVapPress(T-dT,flag,flag2,P,ppmv);
val_b = SatVapPress(T+dT,flag,flag2,P,ppmv);

end

gradient = (val_b - val_a) / (2*dT);
