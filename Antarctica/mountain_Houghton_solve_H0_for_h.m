function [H0]=mountain_Houghton_solve_H0_for_h(Hc,U0,gd)
%finds the critical inversion height for a given mountain height, Hc, U0 and gd
%solving function y=mountain_Houghton_fun_H0_for_h(H0,Hc,U0,gd)
%choose a range for thi in which a solution lies

x=Hc;
sign_fun=sign(mountain_Houghton_fun_H0_for_h(x,Hc,U0,gd));
sign_fun_old=sign_fun;


while sign_fun==sign_fun_old
    x_start=x;  %update the starting location so that fzero doesn't have to do as much work
    x=x+10;
    sign_fun=sign(mountain_Houghton_fun_H0_for_h(x,Hc,U0,gd));
end
    
range=[x_start x];

H0=fzero(@mountain_Houghton_fun_H0_for_h,range,[],Hc,U0,gd);
%function looks like:- y=mountain_Houghton_fun_uAhA(hA,Hc,H0,U0,gd);

