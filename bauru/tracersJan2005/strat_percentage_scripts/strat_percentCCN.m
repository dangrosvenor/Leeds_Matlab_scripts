function [p]=strat_percentCCN(fc)
%works out the percentage increase if convection originally supplies fc % of the troposphere flux of moisture and then increases by 38.5 %
% formula: p=100*(fc*1.385 + (1-fc) -1)

p=100*(fc/100*1.385 + (1-fc/100) -1);