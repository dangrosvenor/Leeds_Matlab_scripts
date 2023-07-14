function y=mountain_smith_func(del,h,L,H0) %

y=del*cos(L*(H0+del-h))-h;   %=0 (eqn 2.13 of Smith (1985)