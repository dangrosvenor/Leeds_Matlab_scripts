function y=mountain_Houghton_fun_uAhA(hA,Hc,H0,U0,gd)

cl=U0 - sqrt(gd*hA/H0*(hA+H0)/2); %eqn (3.8 of Houghton)
uA=(cl*(hA-H0)+H0*U0)/hA;

F0=uA/sqrt(gd*hA);
h_crit = hA .* ( 1 + 0.5*F0.^2 - 1.5*F0.^(2/3) );

y=Hc-h_crit; %h_crit should be consistent with the actual mountain height Hc (i.e. Hc=h_crit, or Hc-H_crit=0)
