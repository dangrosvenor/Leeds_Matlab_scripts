function [Ac,tau,T,SW_out] = calc_SW_surf(f0,W0,N0,SW_in,trans,cf_min,i_Liu)
% Nd in cm3
% W in g/m2
% SW_in - W/m2
% N.B. making assumptions about cloud temperature and pressure in
% albedo_cloudy_func_Seinfeld function
% trans is the estimated transmission ratio of the cloud-free atmosphere (i.e.
% SW_down_surf / SW_down_TOA )

if ~exist('cf_min')
    cf_min = 0.01;
end

if ~exist('i_Liu')
    i_Liu=0;
end

N0 = N0*1e6; %convert to per m3
W0 = W0/1e3; %convert to kg/m2


[Ac,tau] = albedo_cloudy_func_Seinfeld(W0,N0,i_Liu);

Ac(f0<cf_min)=0; %When there is no cloud the cloud albedo will be NaN, so set cloud albedo to zero to return the clear-sky value
%A(A<A_clear) = A_clear; %also don't allow below A_clear (the surface albedo)
        %(for cases when are providing columsn with CF=1 for each one).

%Overall cloud transmission (assume that cloud transmission is 1 minus the
%albedo). And hence when have no cloud the transmission is one. Also, don't need the surface albedo for downwelling SW at surface
T = (1-Ac).*f0 + (1-f0); 
% For when have input of lots of columns :-


SW_out = SW_in .* T .*trans;

%Alternative (equivalent) calculation :-

% SW_surf_cloudy = trans.* SW_in - trans.*SW_in.*Ac;
% SW_surf_clear = SW_in*trans;
% SW_out = f0.*SW_surf_cloudy + (1-f0).*SW_surf_clear; 




