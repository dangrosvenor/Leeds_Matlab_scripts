function y=mountain_Houghton_fun_U0_for_h_H0(U0,Hc,H0,gd)

F0=U0/sqrt(gd*H0);
h_crit = H0 .* ( 1 + 0.5*F0.^2 - 1.5*F0.^(2/3) );

y=Hc-h_crit; %h_crit should be consistent with the actual mountain height Hc (i.e. Hc=h_crit, or Hc-H_crit=0)
