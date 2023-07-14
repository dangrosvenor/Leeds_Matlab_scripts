function y=mountain_Houghton_fun_uBhB(hB,uA,hA,gd)

K3=uA^2/2/gd + hA;
K4=uA*hA;

uB=K4/hB;

y=uB^2/2/gd + hB - K3;
