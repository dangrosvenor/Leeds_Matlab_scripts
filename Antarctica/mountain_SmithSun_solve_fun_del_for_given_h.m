function [del,h]=mountain_SmithSun_solve_fun_del_for_given_h(H0,U,gd)

%gd=N^2*(Hb-Ha)

%range=[-8/L 1.5/L]; %range over which the function should change sign (and thus the solution lie)
range=0;

F0=U/sqrt(gd*H0);
h=H0*( 1 + 0.5*F0.^2 - 1.5*F0.^(2/3) );

del=fzero(@mountain_SmithSun_fun_del_for_given_h,range,[],h,H0,U,gd);  %solve for a given h

%'done mountain wave solve'







