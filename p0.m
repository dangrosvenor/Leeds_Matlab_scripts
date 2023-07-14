function x=p0(T);

th=317;
A=2.53e11;
k=0.286;
B=5.42e3;
T0=th;
rhg=50; 

x=(T/th)^(1/k)*rhg/100*GGEW(T0)-GGEW(T);

%x=A/((T/th)^(1/k)*exp(B/T));

%x=GGEW(T)/(T/th)^(1/k);