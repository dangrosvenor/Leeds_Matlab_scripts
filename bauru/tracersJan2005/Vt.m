function Vt=Vt(af,bf,vis,D,A,k,rhoA,n,g)
%terminal velocity as in Heymsfield paper
Vt=af * (4*g*k/3/rhoA)^bf * vis^(1-2*bf) .* D.^(3*bf-1) * A^((n-1)*bf);