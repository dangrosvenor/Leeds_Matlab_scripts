%top_der.m: returns the derivatives for the symmetric top problem
function ders = top_der(t,w,flag,I,Is,ws,tau0)
%w(1):phi, w(2):phi_dot, w(3):theta, w(4):theta_dot, w(5):psi
%main program produces w(6):psi_dot
ders=[w(2);(Is*ws-2*I*w(2).*cos(w(3))).*w(4)./(I*sin(w(3)));...
       w(4);(tau0-(Is*ws-I*w(2).*cos(w(3))).*w(2)).*sin(w(3))/I;...
       ws-w(2).*cos(w(3));];
