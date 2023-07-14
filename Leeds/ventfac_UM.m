visair=1.44E-5;
Dv=2.26E-5;
vent_1=0.78;
vent_2=0.31; 
mu = 2.5; 
n0 =1246598.69127737381 ;
lam = 8037.63656897062901;
f_x = 0; 
a_x = 130; 
b_x = 0.5;

 

rho=-1.12492679633472403; %error in rho...
rho0=1.21999999999999997;


Sc = visair/Dv;

 T1 = vent_1*(mu+1.)/lam;
    
    T2 = vent_2*Sc.^(1./3.) ...
       * (a_x*(rho/rho0).^(0.5)*rho/visair).^(0.5)
 
    V = 2.*pi*n0/rho*(T1  ...
       + T2*gamma(0.5*b_x + mu + 2.5)/gamma(1.+mu)  ...
       *(1. + 0.5*f_x/lam).^(-(0.5*b_x + mu + 2.5))*lam.^(-0.5*b_x - 1.5))