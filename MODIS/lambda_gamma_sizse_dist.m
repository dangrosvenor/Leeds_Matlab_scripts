function lam = fall_speed_mass_weighted(qR)

%constants for rain
na = 1.1e15;
nb = 0;  %This makes n(D) of the form na * D^alpha * exp(-lambda*D)
c = 523.6;
alpha = 2.5;
d = 3;
rhoa = 1;


lam = ( na*c*gamma(1+alpha+d) ./ (rhoa*qR) ) .^( 1/(1+alpha+d-nb) ) ;



