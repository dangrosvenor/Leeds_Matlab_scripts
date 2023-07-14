function [N,S]=razzak(W,P,T,Na,sig,am)


L=2.5008e6; %latent heat condensation of water
Mw=18e-3; %molecular weight water (kg/mol?)
Ma=28.97e-3; %molecular weight of air
R=8.31; %universal gas constant
Cpa=1005.7;	% specific heat dry air

Ps=SatVapPress(T,'buck2','liq'); %saturation vapour pressure of liquid

t=1;
v=1;
thi=1;
eps=1;
rhoAP=1;
Map=1;

A=2*t*Mw/(1000*R*T);
B=v*thi*eps*Mw*rhoAP/Map/1000; %constants to assign here!

Sm=2/sqrt(B)*(A/3/am).^1.5;
f1=1.5*exp(2.25*log(sig).^2);
f2=1+0.25*log(sig);

alpha=9.81*Mw*L/(Cpa*R*T.^2) - 9.81*Ma/R/T;
gamma=R*T/Ps/Mw + Mw*L^2/(Cpa*P*Ma*T);


Dv=1; %need to define these and possibly constants in Dvd and Kad below.
Kv=1;

%Dvd=Dv/( r/(r+delv) + Dv*(2*pi*Mw/R/T).^0.5/r/alc ); %not sure if need these as are droplet size dependant
%Kad=Ka/( r/(r+delT)+Ka*(2*pi*Mw/R/T).^0.5/r/alT/Cpa );

Dvd=Dv;
Kad=Kv;


G=1000*R*T/(Ps*Dvd*Mw) + L*1000*(L*Mw/R/T - 1)/Kad/T;
G=1./G;

n=(alpha*W/G).^1.5/(2*pi*1000*gamma*Na);
g=2/3*(alpha*W/G).^0.5*A;

S= Sm / ( f1.*(g./n).^1.5 + f2.*(Sm.^2./n).^0.75 ).^0.5;

N=Na*0.5*erfc(  log( f1*(g./n).^1.5 + f2.*(Sm.^2/(n+3.*g)).^0.75 )...
            / (3*sqrt(2)*log(sig))  );



