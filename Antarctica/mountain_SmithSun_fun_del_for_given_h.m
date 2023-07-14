function y=mountain_SmithSun_fun_del_for_given_h(del,h,H0,U,gd)
%NOTE - this is the thin inversion layer approximation

A=1-H0/(H0+del-h);

%y=h + A*H0/(1-A) - del;   %=0  %try different h values until satisfies

%y=del - A*(H0+del) - h*(1-A);

H=H0+del-h;
y=1 - H/H0 + 0.5*U^2/gd/H0*(1-(H0/H)^2);