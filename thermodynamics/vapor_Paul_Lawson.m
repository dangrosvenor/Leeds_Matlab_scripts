function vapor=vapor_Paul_Lawson(tfp) %returns vapour pressure in mb
% INPUT IS IN DEGREES C.  IF GT 0, ASSUMED TO BE DEW POINT.  IF
% LESS THAN 0, ASSUMED TO BE FROST POINT.
% ROUTINE CODES GOFF-GRATCH FORMULA
  tvap=273.16+tfp;
  if(tfp <= 0.) 
% THIS IS ICE SATURATION VAPOR PRESSURE
    if(tvap <= 0) tvap=1E-20; end
    e=-9.09718.*(273.16./tvap-1.)-3.56654.*log10(273.16./tvap)             ...
      +0.876793.*(1.-tvap/273.16);
    vapor=6.1071.*10.^e;
    return
  end
  
% THIS IS WATER SATURATION VAPOR PRESSURE
  if(tvap <= 0) tvap=1E-20; end
    e=-7.90298*(373.16./tvap-1.)+5.02808.*log10(373.16./tvap)             ...
      -1.3816E-7.*(10.^(11.344.*(1.-tvap/373.16))-1.)                  ...
      +8.1328E-3.*(10.^(3.49149.*(1-373.16./tvap))-1);
    vapor=1013.246.*10.^e;
  return