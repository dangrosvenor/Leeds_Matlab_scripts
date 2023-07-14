function derivs = pend2_der( t, x, flag,w0)
% pend2_der: returns the derivatives for the pendulum's full solution
% The function pen2_der describes the equations of motion for a 
% pendulum. The parameter w0, is part of the input
% Entries in the vector of dependent  variables are:
% x(1)-position, x(2)-angular velocity
derivs = [ x(2); -w0^2*sin(x(1))];