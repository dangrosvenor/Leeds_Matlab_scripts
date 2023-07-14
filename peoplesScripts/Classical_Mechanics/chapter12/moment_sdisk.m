%moment_sdisk.m
clear;
f=inline('x.^2.*(1-x.^2).^(0.5)*4/pi'); %define the integrand
J=quad(f,0.0,1.0,1.e-5);                %Simpson quadrature integration 
fprintf('The integral is %4.3f',J);     %print the result to three decimals