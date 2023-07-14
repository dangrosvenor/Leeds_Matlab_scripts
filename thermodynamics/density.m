function [rho]=density(p,T)
%calculates density as ideal gas from pressure and temperature
%function [rho]=density(p,T)
%p in Pa, T in K    

R = 8.314472;

rho=p.*28.97e-3/R./T; %where  28.97e-3 is the molecular weight of air in kg/mol




%%%%%%%%%%%%%%%%%%%%%%%%%%% way did before was like this %%%%%%%%
%M = 28.97*1.67E-27;  %1.67e-27 is the mass of one molecule
% molecular mass of air is mu=28.97e-3 kg/mol. 
% M=mass per molecule (in kg)
% one mole of air weighs 28.97e-3 kg and this is Na=6.022e23 molecules
% so one molecule weighs 28.97e-3/Na. 1e-3/Na = 1.67e-27
% or 1/Na = 1000*1.67e-27 or the weight of one proton in grams

%k = 1.38E-23;
%rho=p*M/k/T;




