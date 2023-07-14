function [Ka,Kad]=thermCond(T,r)

%N= no. droplets per cm^3 and ql= water content g/cm^3

delT=2.16e-5; %cm thermal jump length 
delV=1.096e-5; %cm vapour jump length
alc=1.0;
alT=0.96;
% all above from Abdul-Razzak, 1998

L=2.5008e6; %latent heat condensation of water
Mw=18; %molecular weight water (kg/mol?)
Ma=28.97; %molecular weight of air
R=8.31e7; %universal gas constant as used by Pruppacher and Klett (1980)
Cpa=1005.7;	% specific heat dry air

cal2joul=4.1868;
erg2cal=2.39006e-8; % 1 erg = 1e-7 joules

%r=(ql./N*3/4/pi).^(1/3); %average drop radius in cm

% 1 micron=1e-4 cm
 
%Ka=1.5207e-11.*T.^3 - 4.8574e-8.*T.^2 + 1.0184e-4.*T - 3.9333e-4 ;

% Ka=0.023807 + 7.1128e-5.*(T-273.15); %from Jacobson online lecture notes - same as below with v.slightly different constants

Ka=(5.69+0.017.*(T-273.15))*1e-5; %Pruppacher and Klett, 1977, p.418 - in cal/cm/K/sec
Ka=Ka*cal2joul*100; 

%Dvd=Dv/( r/(r+delv) + Dv*(2*pi*Mw/R/T).^0.5/r/alc ); %not sure if need these as are droplet size dependant

Kad=Ka/( r/(r+delT) + Ka*(2*pi*Mw/R/T).^0.5/r/alT/Cpa ); %Abdul-Razzak

alTpp=0.7;
rho=1e-3;
Cpa=Cpa/1000/cal2joul;
%rho=1;

%R=R/cal2joul;
Ka=5.97e-5;
Kad=Ka/( r/(r+delT) + Ka*(2*pi*Ma/R/T).^0.5/r/alTpp/Cpa/rho ); %Pruppacher

% delVpp=1.3*8e-6;
% Dv=0.3;
% alcpp=0.036
% Kad=Dv/( r/(r+delVpp) + Dv*(2*pi*Mw/R/T).^0.5/r/alcpp );
