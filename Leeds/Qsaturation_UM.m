    function [Qsaturation] = Qsaturation_UM (T, p)    
    %
    % Function to return the saturation mr over water
    % Based on tetans formular
    % QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
    %

    % Temperature in Kelvin
    % Pressure in mb
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
    
    Qsaturation = 999 * ones(size(p));
    
    i = find(T > qsa3 & p .* exp(qsa2 .* (T - tk0c) ./ (T - qsa3) ) > qsa4);
    if length(i)>0
        Qsaturation(i) = qsa1 ./ (p(i) .* exp(qsa2 .* (T(i) - tk0c) ./ (T(i) - qsa3) ) - qsa4);
    end
        
%Can still get zeros here when T approaches qsa3 since if there is a small difference then the expoential amplifies it
% to make the term in the exponential infinite (note that qsa2 and T-tk0c
% are negative) and thus makind Qsaturation zero. Need to make the if
% statement T > qsa3 + small_value. Actually, though, the T-tk0c term also
% makes a big number that causes it to blow up.... So need to limit the
% temperature really.
    
%    i = find(~ (T > qsa3 & p .* exp(qsa2 .* (T - tk0c) ./ (T - qsa3) ) > qsa4) );
%    Qsaturation(i) = 999.;