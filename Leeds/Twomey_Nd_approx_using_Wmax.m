function [Nd,coeff] = Twomey_Nd_approx_using_Wmax(Naero,k,Wmax,w)

%Derivation of % Nd = 0.5 * Naero * w^0.25 :-
% Nd = fA * Naero * fB(k,w)
%k=0.4; Wmax=50;
% fA = 1 / A^(k/(K+2)); 
A =  1.63e-3 * (Wmax*100).^1.5 ./ (k*beta(1.5,k/2));
fA = 1 / A^(k/(k+2)); 
% fB = ( 1.63e-3 * (W*100).^1.5 ./ (k*beta(1.5,k/2)) ) .^(k/(k+2)) 
%    = fC * W.^(1.5*k/(k+2)); % = fC * W.^B
B = (1.5*k/(k+2));
fC = ( 1.63e-3 * (100).^1.5 ./ (k*beta(1.5,k/2)) ) .^(k/(k+2));
coeff = fA.*fC;
% can now calc fB = fC * W.^B and fA
% and use in Nd = fA * Naero * fC * W.^B
% With Wmax=16 m/s and k=0.4 this gives :-
% fA = 0.51, fC = 0.98 and B = 0.25
% fA * fC = 0.5

% Approximation
%Naero=1500;
Nd = coeff.*Naero.*w.^B;