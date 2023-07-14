function derivs = molec_der( t, x, flag)
% molec_der: returns the derivatives for the two atom molecule
% model's full solution
% Entries in the vector of dependent variables are:
% x(1)-position, x(2)-velocity
derivs = [ x(2); 81*pi^2*(3./(x(1).^4)-2./(x(1).^3))/8];