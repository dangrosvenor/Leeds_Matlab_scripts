function qs=qsLES(P,T)

epsilon=0.622;

tk0c = 273.15;      % Temperature of freezing in Kelvin
     qsa1 = 3.8;         % Top in equation to calculate qsat
     qsa2 = -17.2693882; % Constant in qsat equation
     qsa3 = 35.86;       % Constant in qsat equation
     qsa4 = 6.109;       % Constant in qsat equation
     
     
qs = qsa1./(P./100.*exp(qsa2.*(T - tk0c)./(T - qsa3)) - qsa4);

