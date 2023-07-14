function L = LWC_from_h_func(h,T_cb,P_cb)
%returns L in g/m3 for an adiabatic ascent from cloud base temp of T_cb and
%pressure = P_cb over a distance of h (cloud thickness = h)
G = 9.81;
R = 8.314472;

%keep moving upwards from cloud base (in temperature space) until reach the desired h

dT = -0.01; %K
T = T_cb;
P = P_cb;
z = 0;
L = 0;

go=1;
while z<h
    [moist_ad_1,gamm] = moist_ad_lapse_rate(T,P);
%moist_ad_2 is the moist adiabatic temperature without the
%C-C approximation. Also using temperature varation of L.
%moist_ad_1 uses the C-C relationship (with T variation of L)


     rho=P.*28.97e-3/R./T; %where  28.97e-3 is the molecular weight of air in kg/mol
     %rho=p*M/k/T;
     F=-rho*G; %= dp/dz
     [dqldz]=adlwcgm2_just_robs(T,P); %dql/dz in g/m3/m
     
     
     dz =  -1./gamm .* dT;
     z = z + dz;
     T = T +dT;     
     P = P + F.*dz;
     L = L + dqldz.*dz;
                  
end
