function [F]=hydrostatic2(p,z,parr,Tarr)
%solves hydrostatic balance equation in dz/dp form for starting pressure p and height z (m) 
%and temperature profile Tarr (in K as function of pressure, parr)
%same as hydrostatic3 except have made p=-p and z=-z for the c_mex solver
%(needs PSPAN in increasing order)
p=-p;
z=-z;

%M = 28.97*1.67E-27;  %
R = 8.314472;
%k = 1.38E-23;
G = 9.81;

[parr,I]=unique(parr); %remove non-unique pressure values
Tarr=Tarr(I);

[parr,I] = sort(parr); %sort into order
Tarr=Tarr(I);

T=interp1(parr,Tarr,p); %interpolate for temperature for given pressure

rho=p*28.97e-3/R/T; %where  28.97e-3 is the molecular weight of air in kg/mol
%rho=p*M/k/T;
F=-rho*G;
F=1/F'; %dz/dp
'';


