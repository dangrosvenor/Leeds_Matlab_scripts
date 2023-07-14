function y=Smith_dual_layer_solve_del_fun(delA,Ha,Hb,hhat)

delB = - sqrt( delA.^2 + ( (delA - hhat) / (Ha + delA - hhat) ).^2 ); %eqn (23) - try as negative root of delB^2

y = delB.*cos(Hb-Ha+delB-delA) - delA; %eqn (24) of Smith and Sun