function [M] = calc_moment_gamma(p,N,lam,mu)
%function [M] = calc_moment_gamma(p,N,lam,mu)
%Calculates the pth moment of a gamma distribution given number conc,
%lambda and mu. See calc_lambda_n0_gamma_RUN for lambda calculations.

M = N./lam.^p .* (gamma(1+mu+p)) ./ gamma(1+mu); %I guess whether N is in per kg or per m3 will be determined by whether the moment of the n(D) distribution
    %in terms of per kg, or per m3 is required.

