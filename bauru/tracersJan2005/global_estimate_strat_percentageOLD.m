%global estimate of moistening based on the percentage of stratospheric flux of air

%area of globe from -L to +L latitude
L=20; %latitude for tropics (degrees)
R=6400e3; %radius of the Earth (m)

A=4*pi*R^2*sin(pi*L/180); %( integral from -L to +L of dA=2*pi*R'*R*dthi
                   % where R' is the radius of circle at lat thi, R'=Rcos(thi)
                   % dA=circumference of this cicle * arc length R*dthi )
                   
%volume of air entering the strat per year at uplift rate of W m per month
W=1e3; %mean stratospheric uplift (m per month)
V=A*12*W; %volume in one year (12 months) (m^3)

%vapour entering the stratosphere in one year
rho=0.2; %mean density of air at tropopause (kg/m^3)
qvap=4; %mean mixing ratio of air at trop (ppmv)
f=1e6*28.97/18;
Mvap = V * rho * qvap/f %mass of vapour into strat in one year (kg)