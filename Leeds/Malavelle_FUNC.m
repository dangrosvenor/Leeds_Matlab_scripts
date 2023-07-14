function [R,w_fac] = Malavelle_FUNC(dx,Zml)

E1 = 2.59;
E2 = 1.34;
a = 7.95;
b = 8.00;
c = 1.05;

rx=dx./Zml;
R = 1 - (rx.^E1 + a*rx.^E2)./(rx.^E1 + b*rx.^E2 + c); %called sigma* in Malvelle (R in Hamish's draft paper).

w_fac = 2/sqrt(R);
