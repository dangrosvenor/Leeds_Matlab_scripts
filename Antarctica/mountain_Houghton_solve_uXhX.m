function [hA,cl,uA,F0,h_H0,cl_c0]=mountain_Houghton_solve_uXhX(Hc,H0,U0,gd,hB,uB)

%choose a range for thi in which a solution lies

x_start=hB;
x=x_start;
sign_fun=sign(mountain_Houghton_fun_uXhX(x_start,Hc,H0,U0,gd,hB,uB));
sign_fun_old=sign_fun;

while sign_fun==sign_fun_old
    x=x+100;
    sign_fun=sign(mountain_Houghton_fun_uXhX(x,Hc,H0,U0,gd,hB,uB));
end
    
range=[x_start x];

hA=fzero(@mountain_Houghton_fun_uXhX,range,[],Hc,H0,U0,gd,hB,uB);
%function looks like:- y=mountain_Houghton_fun_uXhX(hA,Hc,H0,U0,gd,hB,uB)

cl=U0 - sqrt(gd*hA/H0*(hA+H0)/2); %eqn (3.8 of Houghton)
uA=(cl*(hA-H0)+H0*U0)/hA;
F0=U0/sqrt(gd*H0);
h_H0=Hc/H0;
cl_c0 = cl / sqrt(gd*H0);