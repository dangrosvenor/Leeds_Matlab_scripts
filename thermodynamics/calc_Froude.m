function [F,Hbar]=calc_Froude(U,N,H)
%function [F,Hbar]=calc_Froude(U,N,H)
%Calculates the upstream Froude number, F. Hbar is the dimensionless mountain
%height (=1/F). A high Hbar (>~1.2) suggests blocking in the lower levels.
%See calc_N for function to calculate N (Brunt Vaisala)

F = U ./ (N.*H);
Hbar = 1./F;