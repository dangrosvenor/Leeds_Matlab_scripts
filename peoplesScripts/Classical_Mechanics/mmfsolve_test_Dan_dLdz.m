function y = mmfsolve_test_Dan_dLdz(x)
%x is a vector with values for two different variables (e.g. x(1)=a,
%x(2)=b)

%x(1) = P, x(2) = T, x(3) = z

%Think doesn't need to be a matrix (one value called at a time for each variable)
%These are the funcitons that we want to be zero (check)
y(1) = adlwcgm2_just_robs(x(2),x(1));  %dL/dz  (i.e. dy1/dt)

[ga1,ga2] = moist_ad_lapse_rate(x(2),x(1));
y(2) = -ga2; %dT/dz ( = dy2/dt)

rho = density(x(1),x(2)); G=9.81;
y(3) = -rho.*G; %dP/dz (=dy3/dt)

y=y';

