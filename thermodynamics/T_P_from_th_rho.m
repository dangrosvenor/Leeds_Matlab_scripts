function [T,P]=T_P_from_th_rho(th,rho)
%calculates temperautre (T) and pressure (P) from the potential temperature
%(th) and density (rho) - all in SI units.
%p in Pa, T in K    

R = 8.314472;
C = 28.97e-3/R;

%rho=p.*28.97e-3/R./T; %where  28.97e-3 is the molecular weight of air in kg/mol
%th = T.*( (1000e2./P).^0.286 );
%Can rearrange the above to get :-

k=0.286;
T = (th .* (rho/(1000e2*C)).^k ) .^ (1/(1-k));

P = rho.*T ./ C;