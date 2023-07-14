function f = solve_N_exact_cw_func(H,T_cb,P_cb,tau,re)

%do the integral I
I=quad(@(h)LWC_from_h_ODE(h,T_cb,P_cb,2/3),0,H);

if H==0
    L=0;
    tau2 = 0; %make this part zero, so get a negative f value (otherwise is 0/0=NaN)
else
    [z,y]=ode45(@ode_test_multi_dLdz,[0 H],[0 T_cb P_cb]);
    L = y(end,1);
    rhoW=1000;
    A = 3.*L./(rhoW.*re);
    tau2 = 0.5.*A.*I./(L.^(2/3));
end

rhoW=1000;
A = 3.*L./(rhoW.*re);
f = tau2 - tau;