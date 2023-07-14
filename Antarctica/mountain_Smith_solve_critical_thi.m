function [thi]=mountain_Smith_solve_critical_thi(H0,h,L)

x=0;
sign_fun=sign(mountain_Smith_critical_thi(x,H0,h,L));
sign_fun_old=sign_fun;

while sign_fun==sign_fun_old
    x_start=x;  %update the starting location so that fzero doesn't have to do as much work
    x=x-1/L;
    sign_fun=sign(mountain_Smith_critical_thi(x,H0,h,L));
end
    
range=[x x_start];

del=fzero(@mountain_Smith_critical_thi,range,[],H0,h,L);

thi=H0+del-h;
%function looks like:- y=mountain_Smith_critical_thi(del,H0,h,L)

