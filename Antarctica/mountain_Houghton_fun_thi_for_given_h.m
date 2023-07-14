function y=mountain_Houghton_fun_thi_for_given_h(thi,hx,uA,hA,gd)
%function y=mountain_Houghton_fun_thi_for_given_h(thi,hx,uA,hA,gd)
%finds fluid thickness, thi, for given height of mountain at x, hx,
%and modified upstream velocity and height, uA and hA. And reduced
%gravity, gd.

K3=uA^2/2/gd + hA;
K4=uA*hA;

%y=h + A*H0/(1-A) - del;   %=0  %try different h values until satisfies

%y=del - A*(H0+del) - h*(1-A);

y=K4^2/2/gd/thi^2 + thi + hx - K3; %=0 for solution
% from (3.10) with uc substituted with u=K4/thi from (3.11)