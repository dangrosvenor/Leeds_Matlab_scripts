function [x]=thermAir(P,T,r,choice)

%r=droplet radius in m for kinematic correction
%T=temperautre in K
%P=pressure in Pa
%choice is a string = either 'Ka' or 'Dv' for thermal conductivity or diffusivity of air

P=P/100; %convert P to mb
r=r*100; %convert to cm

Mw=18; %molecular weight water (g/mol)
Ma=28.97; %molecular weight of air
R=8.31e7; %universal gas constant as used by Pruppacher and Klett (1980) i.e. erg/mol/K
Cpa=1005.7;	% specific heat dry air (J/K/kg)

cal2joul=4.1868;


switch choice

case 'Ka'
Ka=(5.69+0.017.*(T-273.15))*1e-5; %Pruppacher and Klett, 1980, p.418 - in cal/cm/K/sec

delT=2.16e-5;
alTpp=0.7;
rho=1e-3; %air density  in g/cm^3
Cpa=Cpa/1000/cal2joul; %heat capacity in cal/g/K required for Pruppacher and Klett

x=Ka/( r/(r+delT) + Ka*(2*pi*Ma/R/T).^0.5/r/alTpp/Cpa/rho ); %Pruppacher
x=x*cal2joul*100; %convert to W/m/K from cal/cm/sec/K

case 'Dv'
delVpp=1.3*8e-6;
Dv=0.211*(T./273.15).^1.94.*(1013.25./P); %in cm^2/sec

alcpp=0.036;
x=Dv/( r/(r+delVpp) + Dv*(2*pi*Mw/R/T).^0.5/r/alcpp ); %correction factor for kinetic effects of droplets of radius r
x=x*1e-4; %convert to m^2/s

end