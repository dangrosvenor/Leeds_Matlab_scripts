qL=2.5e-3;
qL=0;
qv=7e-3;
T_in = 277;
%T_in = 280;

pressure = 879e2; %Pa



Lv = 0.2501e7;   % Latent heat of vapourization
Cp = 1005.; %specific heat capacity
lcrcp = Lv/Cp;



clear qv_save qL_save T_save
    qv_save(1) = qv;
    qL_save(1) = qL;
    T_save(1) = T_in;


for i=1:20
    [cloud_mass_new,T] = Cloud_scheme_test2(qL,qv,T_in);
    dqL = max(cloud_mass_new - qL , -qL);
    qv = qv-dqL;
    qL = qL + dqL;
    T_in = T_in + lcrcp*dqL;  %heat up according to dqL

    qv_save(i+1) = qv;
    qL_save(i+1) = qL;
    T_save(i+1) = T_in;
    
end

qL_save*1e3

T_l = T_in - lcrcp*qL;

qs=Qsaturation_UM(T_in, pressure/100.);
qsL=Qsaturation_UM(T_l, pressure/100.);

RH_tot = (qv+qL)/qs
RH_totL = (qv+qL)/qsL
