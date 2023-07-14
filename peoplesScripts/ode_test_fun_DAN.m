function f=ode_test_fun_DAN(y,t)
%need to return f=dy/dt for vector y and scalar t
% test with dy/dt = y + exp(-t);

f = y + exp(-t);
