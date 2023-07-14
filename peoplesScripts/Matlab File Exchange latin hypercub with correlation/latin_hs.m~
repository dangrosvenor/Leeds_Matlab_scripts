function s=latin_hs(xmean,xsd,nsample,nvar)
% s=latin_hs(xmean,xsd,nsample,nvar)
% LHS from normal distribution, no correlation
% method of Stein
% Stein, M. 1987. Large Sample Properties of Simulations Using Latin Hypercube Sampling. 
%                 Technometrics 29:143-151
% Input:
%   xmean   :  mean of data (1,nvar)
%   xsd     : std.dev of data (1,nvar)
%   nsample : no. of samples
%   nvar    : no. of variables
% Output:
%   s       : random sample (nsample,nvar)
%
% Uses Peter Acklam inverse normal CDF
%
%   Budiman (2003)
% References:
% Iman, R. L., and W. J. Conover. 1980. Small Sample Sensitivity Analysis Techniques for Computer Models, 
% with an Application to Risk Assessment.Communications in Statistics: Theory and Methods A9: 1749-1874
% McKay, M. D., W. J. Conover and R. J. Beckman. 1979.A Comparison of Three Methods for Selecting Values
% of Input Variables in the Analysis of Output from a Computer Code. Technometrics 21: 239-245
%
ran=rand(nsample,nvar);  %random numbers between 0 and 1 of size [nsample nvar]
s=zeros(nsample,nvar);
% method of Stein
for j=1: nvar
   idx=randperm(nsample); %RANDPERM(n) is a random permutation of the integers from 1 to n.
                          %For example, RANDPERM(6) might be [2 4 5 6 1 3].
   P=(idx'-ran(:,j))/nsample;       % probability of the cdf
   s(:,j) = xmean(j) + ltqnorm(P).* xsd(j); % this can be replaced by any inverse distribution function
end

%xsd is the standard deviation, size [1 nvar]

% LTQNORM Lower tail quantile for standard normal distribution.
%  
%     Z = LTQNORM(P) returns the lower tail quantile for the standard normal
%     distribution function.  I.e., it returns the Z satisfying Pr{X < Z} = P,
%     where X has a standard normal distribution.
%  
%     LTQNORM(P) is the same as SQRT(2) * ERFINV(2*P-1), but the former returns a
%     more accurate value when P is close to zero.



