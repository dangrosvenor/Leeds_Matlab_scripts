function P = precip_rate_Wood_2008(W,N)
%function P = precip_rate_Wood_2008(W,N)
% Returns precipitation rate in mm/hour given the LWP (W) in g/m2
% and N in cm^-3 using the formula
% for precip rate given in Wood, JGR (2008). (Title :- Open celluar...)
% Also see re_precip_rate.m function to get the reff for a given tau and
% precip rate (i.e. inverting the relationship here).

W(W<0)=0;
N(N<5)=5;
P = (1/24) * 0.37 * (W(:)./N(:)).^1.75;  %converting here to mm/hr from mm/day as originally in Wood paper.
P = reshape(P,size(W));

%Bennartz (2007) gives an alternative formula of P propto W^3/2 / N
%For comparison - CPT precip rates ranged from 0.01 to 0.03 mm/hr moving
%from east to west (mean values, though).

