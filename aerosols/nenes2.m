function [I]=nenes2(smax,P,T,W,N,sc)
%finds smax as for Nenes and Seinfeld parameterisation
%sc= supersaturation vector in % above 100%.


sc=sc/100; %convert to e/es form


if smax<=0;smax=1e-10;end

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

A=4*Mw*tensW/R/T/rhoW;


%A=ones([1 length(sc)])*A;

Ps=SatVapPress(T,'buck2','liq'); %saturation vapour pressure of liquid in Pa
Dv=thermAir(P,T,1e2,'Dv');
Ka=thermAir(P,T,1e2,'Ka');

G=rhoW*R.*T./(Ps*Dv*Mw) + L*rhoW./(L*Mw/R./T - 1)/Ka./T;

G=rhoW*R.*T./(Ps*Dv*Mw) + L*rhoW.*(L*Mw/R./T - 1)/Ka./T;
G=4./G;

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
goA=1;
goB=1;

del=smax^4 - 16*A^2*alpha*W/9/G;

if del>=0
    spart=0.5*(1 + ( 1 - 16*A*A*alpha*W/(9*G*smax^4) )^0.5 );
    spart=sqrt(spart)*smax;
else
    spart=smax * min([ 2e7*A/3*smax^-0.3824 1.0]);
end


%s needs to be less than smax
    %for asin(s/smax) and sqrt(smax^2-s^2) for first half of integration
    
    
[mindiff imax]=min(abs(sc-smax));
if length(imax)>1; imax=imax(1); end
if sc(imax)<smax; imax=imax+1; end %in case closest is on lower side (then would miss last point)
if imax>length(sc); smax=sc(end); imax=length(sc); end %smax can't go above sc(end)
if imax==1; goA=0; goB=0; Ia=0; Ib=0; end %smax smaller than min sc so abort and return zero integrals (if was originally between bin 1&2 imax now =2)

[mindiff ipart]=min(abs(sc-spart)); %find ipart as index of sat nearest to spart
if length(ipart)>1; ipart=ipart(1); end
if ipart==imax; ipart=ipart-1; end %if nearest is on the upper side and =imax (sc(imax)>smax)
if ipart==1; ipart=2; Ia=0; goA=0; end %if ipart is on lower side bwetween 1 and 2 (can't be 1) aborts Ia

%if ipart>imax; goB=0; Ib=0; end

if goA==1
for j=1:ipart-1 %added one to s index since in paper s runs from 0 to imax
    
    x=sc(j+1);
    Iaa2=x/2 .* sqrt(smax^2 - x.^2) + smax.^2/2 .* asin(x./smax);
    x=sc(j);
    Iaa1=x/2 .* sqrt(smax^2 - x.^2) + smax.^2/2 .* asin(x./smax);
    
    suma(j)= N(j)/(sc(j+1)-sc(j)) * (Iaa2-Iaa1);
    
end
end

if exist('suma')==0; suma=0; end
Ia=sqrt(G/alpha/W) * sum(suma);


A=ones([1 length(sc)])*A;


if goB==1
for j=ipart-1:imax-2
    
	sumb= N(j)*A(j)/(sc(j+1)-sc(j)) * log(sc(j+1)./sc(j)) + (N(imax-1)*A(imax-1)/( sc(imax)-sc(imax-1) )) * log(smax/sc(imax-1)); %changed imax to imax-1 for N
end

if exist('sumb')==0; sumb=0; end
Ib=2/3 * sum(sumb);

end

I=pi/2*gamma*rhoW*G*smax/alpha/W * ( Ia + Ib ) - 1; %this needs to be zero - guess smax

spart=spart/smax;
%smax