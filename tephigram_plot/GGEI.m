function ei=GGEI(T)


% goff gratch formulae
Y=273.16./T;
ei=6.1071.*(10.^(-9.09718.*(Y-1.)...
     		   -3.56654.*log10(Y)...
     	   +0.876793.*(1.-(1./Y))));

ei=ei.*100;