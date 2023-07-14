function [I,spart,smax]=spart(P,T,W,N,sc,smax)
%finds smax as for Nenes and Seinfeld parameterisation
%sc= supersaturation vector in % above 100%.


sc=sc/100; %convert to e/es form

L=2.5008e6; %latent heat condensation of water
Mw=18.016e-3; %molecular weight water (kg/mol?)
Ma=28.9644e-3; %molecular weight of air
R=8.3144; %universal gas constant
Cpa=1005.7;	% specific heat dry air
rhoW=997; %density of water


mu=3; %no. ions of salt that dissociate
Ms=132.1e-3; %molecular weight of salt kg/mol
rhoS=1.769e3; %density of salt kg/m^3 %ammonium sulphate=1.77e3
tensW=0.07197; %J/m^2 %surface tension water
                %variation with temperature?
%tensW=115e-3;

%tensW=pi/2*tensW;



%A=ones([1 length(sc)])*A;

Ps=SatVapPress(T,'buck2','liq'); %saturation vapour pressure of liquid in Pa
Dv=thermAir(P,T,1e2,'Dv');
Ka=thermAir(P,T,1e2,'Ka');

G=rhoW*R.*T./(Ps*Dv*Mw) + L*rhoW./(L*Mw/R./T - 1)/Ka./T;

G=rhoW*R.*T./(Ps*Dv*Mw) + L*rhoW.*(L*Mw/R./T - 1)/Ka./T;
G=4./G;

A=4*Mw*tensW/R/T/rhoW;


alpha=9.81*Mw*L/(Cpa*R*T.^2) - 9.81*Ma/R/T;
gamma=P*Ma/Ps/Mw + Mw*L^2/(Cpa*R.*T.^2);

%alpha=alpha/100; %convert to cm^-1?????

% I=pi/2*gamma*rhoW*G*smax/alpha/W * ( Ia + Ib ) - 1; %this needs to be zero - guess smax
% 
% 
% 
% Iaa=x/2 .* sqrt(smax^2 - x.^2) + smax.^2/2 .* arcsin(x./smax); %do for x=all scj j=0:ipart
% 
% suma= Na ./ (sca2-sca1) .* ( Iaa(2:end) - Iaa(1:end-1) )
% 
% Ia=sqrt(G/aplha/W) * sum(suma);
% 
% 
% sumb= N.*A./(scb2-scb1) .* log(scb2./scb1) + (Nimax*Aimx/( scb2(end)-scb1(end) )) * log(smax/scb1(end));
% 
% Ib=2/3 * sum(sumb);


del=smax.^4 - 16*A^2*alpha*W/9/G;

a=find(del>=0);
    spart(a)=0.5*(1 + ( 1 - 16*A*A*alpha*W./(9*G.*smax(a).^4) ).^0.5 );
    spart(a)=sqrt(spart(a)).*smax(a);
    
b=find(del<0);
ss=2e7*A/3.*smax(b).^-0.3824;
bb=find(ss>1.0);
ss(bb)=1.0;

    spart(b)=smax(b) .* ss;
   

spart=spart./smax;