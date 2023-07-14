function s=latin_hs_owndist_Dan(nsample,nvar,x_pdf,f_pdf)
% function s=latin_hs_owndist_Dan(nsample,nvar,x_pdf,f_pdf)
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
% Dan - added the option of supplyin own PDF instead of a gaussian with
% given std deviation. x_pdf{nvar}(1:npdf_vals) is an array of bin mid-points
% and f_pdf{nvar}(1:npdf_vals) are the frequencies.


ran=rand(nsample,nvar);  %random numbers between 0 and 1 of size [nsample nvar]
s=zeros(nsample,nvar);
% method of Stein
for j=1: nvar
   idx=randperm(nsample); %RANDPERM(n) is a random permutation of the integers from 1 to n.
                          %For example, RANDPERM(6) might be [2 4 5 6 1 3].
   P=(idx'-ran(:,j))/nsample;       % probability of the cdf

   %the distribution for P(x) at this stage has mean=0 and std_dev=1 (a
   %standard normal distribution).
% s(:,j) = xmean(j) + ltqnorm(P).* xsd(j); % this can be replaced by any inverse distribution function   
     %this line converts from a standard normal distribution (i.e. mean=0
     %and std=1) to a normal dist with a given mean and std - since ltqnorm
     %below is for a standard normal distribution. The variable transform
     %to covert from a non-standard normal to a standard normal is
     %is x2 = (x-mean_x)/std_x. This does not change the factor in front of the 
     %intergral over this distribution. So above it is just converting back
     % i.e. x = std_x*x2 + mean_x
     
     %To add our own distribution we just need to supply the x-value, i.e.:-          
        s(:,j) = PDFdistZ(x_pdf{j},f_pdf{j},P);

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



