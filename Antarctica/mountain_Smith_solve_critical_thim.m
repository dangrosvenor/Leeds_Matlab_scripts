function [U0]=mountain_Houghton_solve_U0_for_h_H0(Hc,H0,gd)

%solving function function y=mountain_Houghton_fun_U0_for_h_H0(U0,Hc,H0,gd)
%choose a range for thi in which a solution lies

x=0;
sign_fun=sign(mountain_Houghton_fun_U0_for_h_H0(x,Hc,H0,gd));
sign_fun_old=sign_fun;


while sign_fun==sign_fun_old
    x_start=x;  %update the starting location so that fzero doesn't have to do as much work
    x=x+1;
    sign_fun=sign(mountain_Houghton_fun_U0_for_h_H0(x,Hc,H0,gd));
end
    
range=[x_start x];

U0=fzero(@mountain_Houghton_fun_U0_for_h_H0,range,[],Hc,H0,gd);
%function looks like:- y=mountain_Houghton_fun_uAhA(hA,Hc,H0,U0,gd);

