function B = Nd_equations_B_func_of_tau_W(CTT)
% function B = Nd_equations_N_func_of_tau_W()
% Calculates the constant in front of equation for N as function of tau and
% LWP (=W); N = B * tau^3 / W^(5/2)

if ~exist('CTT')
    CTT=280; %pick a sensible value
end

%set a tau and re
tau_set=20;
re_set = 12e-6;

%This a function of CTT only for this calc (P is assumed in
%MODIS_justN_func
N=MODIS_justN_func(tau_set,re_set,'calc',0,CTT,'N');
W=MODIS_justN_func(tau_set,re_set,'calc',0,CTT,'W');


B = N .* W.^(5/2) ./ tau_set.^3;

