function y=mountain_Houghton_fun_H0_for_h(H0,Hc,U0,gd)
%function y=mountain_Houghton_fun_H0_for_h(H0,Hc,U0,gd)
%calculates the critical H0 for a given mountain height, Hc and U0 and gd

F0=U0/sqrt(gd*H0);
h_crit = H0 .* ( 1 + 0.5*F0.^2 - 1.5*F0.^(2/3) );

y=Hc-h_crit; %h_crit should be consistent with the actual mountain height Hc (i.e. Hc=h_crit, or Hc-H_crit=0)
