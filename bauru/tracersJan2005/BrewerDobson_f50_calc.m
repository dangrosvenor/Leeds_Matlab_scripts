% 1+f50 = 1.25*(1 + f50/1.385)
%       = 1.25 + f50*1.25/1.385
% f50*(1 - 1.25/1.385) = 0.25

p=22.5; %percentage of trend of vapour trying to describe
pi=108; %38.5 for abv trop, 108 % for below %percentage that convective portion increases moistening by

f50 = p/100 / (1 - (1+p/100)/(1+pi/100) )

per = 100 * f50/(1 + f50)  %percentage that convection then represents of total flux