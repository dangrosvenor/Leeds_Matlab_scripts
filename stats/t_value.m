function [t,serror,t2,serror2,Neff] = t_value(trend,residuals,x,y)
%function [t] = t_value(trend,residuals,x)
% Calculates the t-value for the linear fit y(x) = ax + c as a function of the resdiuals
% and the x-values (NOT y-values)!

%Got the following from Chapter 1 of some lecture notes I found -
%Chapter 1 - Data Analysis :-
%<file://C:\Users\eardgro\Documents\logbook_OLD\lecture notes etc\Statistical_significance_of_linear_trend.pdf>
%Is based on Santer et al. - Statistical significance of trends and trend
%differences in layer-average atmospheric temperature time series,
%JGR, 2000

N=length(x);

se_sq = sum(residuals.^2)/(N-2);
me = mean(x); %N.B. - these are the x values (not y)!
sum_sq = sum( (x - me).^2 );
sa_sq = se_sq ./ sum_sq;
serror = sqrt(sa_sq);
t = trend/sqrt(sa_sq);

%rau = autocorrel_Dan(y);
rau = autocorrel_Dan(residuals); %Do auto-corr of residuals not the actual values, following Santer (2000)
r1 = rau(2); %Lag-1 autocorrelation coefficient - N.B., this can be <0, which increases the 
%effective sample size - not sure what to in that case. Perhaps best to use
%the smallest of the sample sizes?
Neff = N * (1-r1)./(1+r1); %Eqn. 6 of Santer (2000).
se_sq2 = sum(residuals.^2)/(Neff-2); %Eqn. 4 of Santer (2000).
sa_sq2 = se_sq2 ./ sum_sq; %Eqn. 3. of Santer (2000).
serror2 = sqrt(sa_sq2);
t2 = trend/sqrt(sa_sq2);

  