function [F]=hydrostatic(z,p,zarr,Tarr)
%solves hydrostatic balance equation for pressure from starting pressure p0 and height z0 (m) and temperature profile T (in K)

M = 28.97*1.67E-27;
k = 1.38E-23;
G = 9.81;
T=interp1(zarr,Tarr,z); %interpolate for temperature for given z

rho=p*M/k/T;
F=-rho*G;
F=F';



