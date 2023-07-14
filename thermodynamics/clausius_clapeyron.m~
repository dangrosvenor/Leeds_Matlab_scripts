%assuming constant pressure
clear e

dT=1e-3;

R=287.04;    
L=2.5e6;  
eps=0.622;

T0=273;
T1=300;

%N = (T1-T0)/dT;
T=T0:dT:T1;
N=length(T);


e(1) = SatVapPress(T(1),'goff','liq');

for i=2:N
    
    e(i) = e(i-1) + L*e(i-1)/(R*T(i-1).^2/eps) * dT;
    
end

%more complicated/accurate formula
e2 = SatVapPress(T,'goff','liq');
dedT2 = diff(e2)/dT;

dedT = L*eps .*e2 /R ./T.^2;
%actually C-C is much closer if we take into account the temperature
%variation of L - use latent_heat_temp_variation.m


%theoretical solution of integral of C-C
e3 = e(1)*exp(-L*eps/R .*(1./T - 1./T(1)));
%i.e. C-C is just an approximation for the saturation vapour pressure -
%likely to be less accurate than the more complicated formulae
%they diverge more at warmer temperatures

