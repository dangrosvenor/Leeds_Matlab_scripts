function [qV,qL,qT,qs,qsL,RHtot,dmass] = read_cld_scheme_output_func(read_file,format_str,cloud_scheme);


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

%index 12 is the first lot of data

qT = dat{13}(:);
qs = dat{14}(:);
qsL = dat{15}(:);
T = dat{16}(:);
TL = dat{17}(:);
qL = dat{18}(:);
RHcrit = dat{19}(:);
docloud = dat{20}(:);

RHtot = qT./qsL;


% do calculation of what the original saturation adjustment scheme would have calculated

if cloud_scheme==1
    qV = qT - qL; %reverse calculate qv from total q
else
    qV = qT;
end


dqsdt = dqwsatdt_fun(qs, T);
qsatfac = 1./(1. + Lv/Cp*dqsdt);

dmass = max ( -qL, (qV - qs) .* qsatfac );
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
         
         