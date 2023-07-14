function y=mountain_Houghton_fun_uXhX(hX,Hc,H0,U0,gd,hB,uB)

cr=uB - sqrt(gd*hX/hB*(hX+hB)/2); %eqn (3.8 of Houghton)
uX=(cr*(hX-hB)+hB*uB)/hX;

y = uX - 2*sqrt(gd*hX) - U0 + 2*sqrt(gd*H0);

