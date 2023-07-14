function ew=GGEW(T)


% goff gratch formulae - returns sat vapour pressure in Pa. T in K
Y=373.16./T;
ew=1013.246.*10.0.^(-7.90298.*(Y-1.0)...
     		    +5.02808.*log10(Y)...
     		    -1.3816E-07.*((10.0.^(11.344.*(1.0-(1.0./Y))))-1.0)...
     		    +8.1328E-03.*((10.0.^(-3.49149.*(Y-1.0)))-1.0));
            
ew=ew.*100;
