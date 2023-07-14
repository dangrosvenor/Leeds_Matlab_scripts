function [k,qV,qL,qT,qs,dqsdt,qsatfac,dt,qL_new,dqL] = read_cld_scheme_output_postcond_func(read_file,format_str,cloud_scheme);


fid=fopen(read_file,'rt'); %is quicker to use fscanf that dlmread
dat = textscan(fid,format_str);
fclose(fid);


%chop off the last line as can be incomplete
nrec = -1;
for i=1:length(dat)
    nrec = max(nrec,length(dat{i}) );
end
for i=1:length(dat)
    dat{i} = dat{i}(1:nrec-1);
end




Lv = 0.2501e7;   % Latent heat of vapourization
Cp = 1005.; %specific heat capacity

dat_col01 = 12; %index for first column of data

icol=0;
k = dat{dat_col01 + icol}(:); icol=icol+1; 
qT = dat{dat_col01 + icol}(:); icol=icol+1; 
qs = dat{dat_col01 + icol}(:); icol=icol+1; 
qL = dat{dat_col01 + icol}(:); icol=icol+1; 
dqL = dat{dat_col01 + icol}(:); icol=icol+1; 
dqsdt = dat{dat_col01 + icol}(:); icol=icol+1; 
qsatfac = dat{dat_col01 + icol}(:); icol=icol+1; 
dt = dat{dat_col01 + icol}(:); icol=icol+1; 
qL_new = dat{dat_col01 + icol}(:); icol=icol+1; 


% do calculation of what the original saturation adjustment scheme would have calculated

if cloud_scheme==1 %if cloud scheme is on then qT as read in is qtotal, otherwise it is qV
    qV = qV - qL; %reverse calculate qv from total q
else
    qV = qT;
    qT = qV + qL;
end


%dqsdt = dqwsatdt_fun(qs, T);
%qsatfac = 1./(1. + Lv/Cp*dqsdt);

%dmass = max ( -qL, (qV - qs) .* qsatfac );
    % -qL is the lower bound of dqL. Otherwise will remove more than have.
    % if dqL is positive then this is not relevant
'';



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
    
    
    function [qsat] = Qsaturation (T, p)    
    !
    ! Function to return the saturation mr over water
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
    !

    ! Temperature in Kelvin
    ! Pressure in mb
    tk0c = 273.15;
    qsa1 = 3.8;
    qsa2 = - 17.2693882;
    qsa3 = 35.86;
    qsa4 = 6.109;
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
    
    i = find(T > qsa3 & p .* exp(qsa2 .* (t - tk0c) ./ (T - qsa3) ) > qsa4);
    Qsaturation(i) = qsa1 / (p(i) .* exp(qsa2 .* (T(i) - tk0c) / (T(i) - qsa3) ) - qsa4);
    
    i = find(~ (T > qsa3 & p .* exp(qsa2 .* (t - tk0c) ./ (T - qsa3) ) > qsa4) );
    Qsaturation(i) = 999.



    
    
    
         
         