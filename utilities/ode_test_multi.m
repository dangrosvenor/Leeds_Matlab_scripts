function F = ode_test_multi(t,y)
%y(1) = L, y(2) = T, y(3) = P, t = z (the stepping variable - known)

%call using [t,y]=ode45(@ode_test_multi,[t0 t1],[y1_0 y2_0]);

F(1,1) = adlwcgm2_just_robs(y(2),y(3));  %dL/dz  (i.e. dy1/dt)

[ga1,ga2] = moist_ad_lapse_rate(y(2),y(3));
F(2,1) = -ga2; %dT/dz ( = dy2/dt)

rho = density(y(3),y(2)); G=9.81;
F(3,1) = -rho.*G; %dP/dz (=dy3/dt)