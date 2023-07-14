function [cr]=mountain_Houghton_solve_cr(U0,Hc,H0,gd)

[hA,cl,uA,Fa,xa,Cl,uB,hB,FB,DB,hX,FX,DX,Cr]=mountain_Houghton_solve_uBhB(Hc,H0,U0,gd);

cr=Cr*sqrt(gd*H0);