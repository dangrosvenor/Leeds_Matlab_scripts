cx(1)=523.6; %rain
cx(2)=52.36; %snow
cx(3)=366.5;   %261.8; - changed to this due to tropical graupel modification
cx(4)=104; %ice
%N.B. this cx value doesn't actually affect the Z calc, since it cancels out

dx=3;

nax(1)=1.1e15;
nax(2)=2.0e27;
nax(3)=5e25;

nbx(1)=0;
nbx(2)=-3.5;
nbx(3)=-4;

alx(1)=2.5;
alx(2)=2.5;
alx(3)=2.5;
alx(4)=0;

fx(1)=0; %fall speed f-values are all zero
fx(2)=0;
fx(3)=0;
fx(4)=0;

bx(1)=0.65; %fall speed b-values 
bx(2)=0.25;
bx(3)=0.734;
bx(4)=0.6635;

ax(1)=362; %fall speed a-values 
ax(2)=4.84;
ax(3)=253;
ax(4)=71.34;


qn=2.4e8;
ql=[0:0.1:5]/1e3;
ql=5e-3;
qr=[0:0.1:5]/1e3;

qlem=0.001*(ql - (pi/6*1000*qn*20e-6^3/rho) );
qmorr=1350*ql.^2.47 * (qn*rho/1e6)^-1.79;

phys1d_25 = cx(1)*nax(1)* gamma(1+alx(1)+dx)/rho ;
phys1d_64 = 1*pi/4 * ax(1) * gamma(3 + alx(1) + bx(1) ) * (1.22/rho)^0.5 ; 

Lamb_R = (phys1d_25./qr).^(1./(1+dx+alx(1)));
%N0_R = nax(1)*Lamb_R^nbx(1);
N0_R = nax(1); %since nbx(1)=0
PRaciTEMP = phys1d_64 * N0_R * (Lamb_R + fx(1)).^-(alx(1)+3+bx(1));
pracw_lem = 1/1 * PRaciTEMP .* ql ;

pracw_morr = 67 * (ql*qr).^1.15;




