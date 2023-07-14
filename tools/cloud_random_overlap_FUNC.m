function [CFran, CFmin, CFmax] = cloud_random_overlap_FUNC(f)
% function [CFran, CFmin, CFmax] = cloud_random_overlap_FUNC(f)
% f is a profile of cloud fractions.
% Returns the random overlap of the profile based on :-
% https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/JD094iD07p09925
% (Tian and Curry, JGR, 1989).
% Also returns the minimum overlap and max overlap values


CFran=1;
for i=1:length(f)
    CFran = CFran - (1-f(i));
end

CFmin = min(sum(f),1);
CFmax = max(f);