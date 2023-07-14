function [dqwsatdt] = dqwsatdt_fun(qsat,T)
    % Function to return the rate of change with temperature
    % of saturation mixing ratio over liquid water.
    %
    % Based on tetans formular
    % QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  %T is temp in K, saturation MR (kg/kg)

    tk0c = 273.15;
    qsa1 = 3.8;
    qsa2 = - 17.2693882;
    qsa3 = 35.86;
    qsa4 = 6.109;
    % Temperature of freezing in Kelvin
    % Top in equation to calculate qsat
    % Constant in qsat equation
    % Constant in qsat equation
    % Constant in qsat equation
    %
    
    dqwsatdt = - qsa2 .* (tk0c - qsa3) .* (1.0 + qsa4 * qsat / qsa1) .* qsat .* (T - qsa3) .^ ( - 2.0);