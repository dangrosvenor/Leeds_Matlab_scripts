%torquef2_der.m: returns the derivatives for the torque free ellipsoid
function derivs = torquef2_der( t, w, flag, gam1, gam2, gam3)
% w(1):w1, w(2):w2, w(3):w3
derivs =[-w(2).*w(3)*gam1;-w(3).*w(1)*gam2;-w(1).*w(2)*gam3;];

