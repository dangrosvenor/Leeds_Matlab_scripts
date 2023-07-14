%finds saturation temp for given RH and ground conditions (assumes dry adiabatic rise)

function y=Ts

y=fzero(@zero,273);


function x=zero(T);

th=317;
A=2.53e11;
k=0.286;
B=5.42e3;
T0=th;
rhg=50;
p=(T/th)^(1/k); %actually =p/p0


x=p*rhg/100*GGEW(T0)-GGEW(T);
