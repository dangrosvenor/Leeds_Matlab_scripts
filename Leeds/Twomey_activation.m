%Twomey 1959 (part 2) paper has a maritime and continental aerosol
%distribuiton (values given under Table 1). 
% Continental : c=2000, k=0.4
% Maritme: c=125, k=1/3

% Number of activated droplets based on N=CS^k where S is the max
% supersat (Eqn. 11).
% Eqn. (10) gives Smax as function of V (updraft in cm/s).

% But need to relate this to a given aerosol amount. Hard to do this
% properly, but Paul's CASIM code assumes that all the aerosol is activated
% by the time we reach a certain updraft Wmax=16 m/s and k=0.4. Then can use Eq. 10
% to calculate the Smax at which all aerosol activates. Call this Sall.
% Then can use Naero=C*Sall^k to calculate C for this distribution.
% Doing this (converting to use m/s for W) the formula for Nd boils down to
% Nd = 0.5 * Naero * w^0.25
% Quick check here to see if it is similar to Twomey using the maritmie and
% continental c and k values and proper relationship (Eqn. 11)

%Derivation of % Nd = 0.5 * Naero * w^0.25 :-
% Nd = fA * Naero * fB(k,w)
k=0.4; Wmax=50;
% fA = 1 / A^(k/(K+2)); 
A =  1.63e-3 * (Wmax*100).^1.5 ./ (k*beta(1.5,k/2));
fA = 1 / A^(k/(k+2)); 
% fB = ( 1.63e-3 * (W*100).^1.5 ./ (k*beta(1.5,k/2)) ) .^(k/(k+2)) 
%    = fC * W.^(1.5*k/(k+2)); % = fC * W.^B
B = (1.5*k/(k+2));
fC = ( 1.63e-3 * (100).^1.5 ./ (k*beta(1.5,k/2)) ) .^(k/(k+2));
fA*fC
% can now calc fB = fC * W.^B and fA
% and use in Nd = fA * Naero * fC * W.^B
% With Wmax=16 m/s and k=0.4 this gives :-
% fA = 0.51, fC = 0.98 and B = 0.25
% fA * fC = 0.5


figure

W = [0.1:0.1:16];

% Continental : 
c=2000; k=0.4;
N = c.^(2./(k+2)) .* (1.63e-3.*(W.*100).^(3/2) ./ (k.*beta(1.5,k/2)) ) .^(k/(k+2));
plot(W,N,'b','linewidth',2);
hold on
leg{1} = 'Twomey continental';

% Maritme: 
c=125; k=1/3;
N = c.^(2./(k+2)) .* (1.63e-3.*(W.*100).^(3/2) ./ (k.*beta(1.5,k/2)) ) .^(k/(k+2));
plot(W,N,'r','linewidth',2);
hold on
leg{2} = 'Twomey maritime';

% Approximation
Naero=1500;
N = fA*fC.*Naero.*W.^B;
plot(W,N,'b--','linewidth',2);
hold on
leg{3} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero)];

% Approximation
Naero=100;
N = fA*fC.*Naero.*W.^B;
plot(W,N,'r--','linewidth',2);
hold on
leg{4} = ['Approx wmax=' num2str(Wmax) ', Naero=' num2str(Naero)];

legend(leg);

