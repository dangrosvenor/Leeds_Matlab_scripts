function [M_over_ND]=Moments_for_wood2005_dist(n,rmean,r0)
%Returns the nth moment of the truncated distribution given in Wood, 2005,
%Drizzle Part 2 divided by ND (i.e. M(n)/ND )
%Then can multiply by Nd to get teh moment, or use it to calculate ND from
%.e.g. the rain rate.
%r0 was set to 20 microns in the Wood paper.
% COmtaock paper found that rmean was fairly constant at 40um, so one approach might be
% to assume this and estimate ND from this function and the rain rate using
% R = 4/3*pi * AT * rhow * M(d+3) where the fall speed relation v(r) =
% AT*r^d has been used.

ivals = 0:floor(n)+5;

for ii=ivals
    MM(ii+1) = Moments_Wood_sum_part(ii,rmean,r0);
end

M_over_ND = interp1(ivals,MM,n);
%M_over_ND = 10.^ (interp1(ivals,log10(MM),n) ); %Or do in log10 space??
%M_over_ND = interp1(ivals,MM,n,'pchip');
%M_over_ND = interp1(ivals,MM,n,'spline');
  %Result is very sensitive to how this interpolation is done...

function [MM] = Moments_Wood_sum_part(n,rmean,r0)

sum=0;
for i=0:n
    sum = sum + (r0./(rmean-r0)).^i ./ (factorial(i));
end
MM = factorial(n) .* (rmean-r0).^n .* sum;


