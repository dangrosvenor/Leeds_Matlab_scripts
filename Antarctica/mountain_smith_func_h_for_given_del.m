function y=mountain_smith_func_h_for_given_del(h,del,L,H0)

y=del*cos(L*(H0+del-h))-h;   %=0  %try different h values until satisfies