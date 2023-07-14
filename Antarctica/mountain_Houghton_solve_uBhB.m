function [hA,cl,uA,F0,h_H0,cl_c0,uB,hB,FB,DB,hX,FX,DX,Cr]=mountain_Houghton_solve_uBhB(Hc,H0,U0,gd)

%choose a range for thi in which a solution lies

sign_fun=sign(mountain_Houghton_fun_uAhA(H0,Hc,H0,U0,gd));
sign_fun_old=sign_fun;
x=H0;
while sign_fun==sign_fun_old
    x_start=x;
    x=x+100;
    sign_fun=sign(mountain_Houghton_fun_uAhA(x,Hc,H0,U0,gd));
end
    
range=[x_start x];

hA=fzero(@mountain_Houghton_fun_uAhA,range,[],Hc,H0,U0,gd);
%function looks like:- y=mountain_Houghton_fun_uAhA(hA,Hc,H0,U0,gd);

cl=U0 - sqrt(gd*hA/H0*(hA+H0)/2); %eqn (3.8 of Houghton)
uA=(cl*(hA-H0)+H0*U0)/hA;
F0=U0/sqrt(gd*H0);
h_H0=Hc/H0;
c0=sqrt(gd*H0);
cl_c0 = cl / c0;

x=hA-100;
sign_fun=sign(mountain_Houghton_fun_uBhB(x,uA,hA,gd));
sign_fun_old=sign_fun;

while sign_fun==sign_fun_old
    x_start=x;
    x=x-100;
    sign_fun=sign(mountain_Houghton_fun_uBhB(x,uA,hA,gd));
end
    
range=[x x_start];
hB=fzero(@mountain_Houghton_fun_uBhB,range,[],uA,hA,gd);
%y=mountain_Houghton_fun_uBhB(hB,uA,hA,gd)
uB=uA*hA/hB;

FB=uB/c0;
DB=hB/H0;


[hX]=mountain_Houghton_solve_uXhX(Hc,H0,U0,gd,hB,uB);
cr=uB - sqrt(gd*hX/hB*(hX+hB)/2); %eqn (3.8 of Houghton)
uX=(cr*(hX-hB)+hB*uB)/hX;

FX=uX/c0;
DX=hX/H0;
Cr=cr/c0;



