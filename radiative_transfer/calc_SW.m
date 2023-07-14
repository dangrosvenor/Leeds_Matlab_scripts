function [Ac,tau,A,SW_out] = calc_SW(fc,Wc,Nc,SW_in,A_surf,trans,cf_min,iprovide_tau,tau,i_Liu)
% Calculates the SW upwelling at TOA based (SW_out) on the cloud properties
% and SW downwelling at cloud top (SW_in) and atmospheric transmissivity.
% fc is the cloud fraction in the range of 0 to 1.
% Nd in cm3
% W in g/m2
% SW_in - downward SW AT CLOUD TOP in W/m2 - can be estimated using SW down at the
% surface for clear skies. Or SW_down TOA * transmissivity.
% N.B. making assumptions about cloud temperature and pressure in
% albedo_cloudy_func_Seinfeld function
% trans is the estimated transmission ratio of the atmosphere (i.e.
% SW_down_surf / SW_down_TOA )
% A_surf is the surface albedo.

if ~exist('cf_min')
    cf_min = 0.01;
end
if ~exist('iprovide_tau')
    iprovide_tau = 0;
end
if ~exist('i_Liu')
    i_Liu=0;
end

Nc = Nc*1e6; %convert to per m3
Wc = Wc/1e3; %convert to kg/m2

if iprovide_tau==0
    [Ac,tau] = albedo_cloudy_func_Seinfeld(Wc,Nc,i_Liu);
else
    Ac = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
end

tau(fc<cf_min)=NaN;
Ac(fc<cf_min)=0; %When there is no cloud the cloud albedo will be NaN, so set cloud albedo to zero to return the clear-sky value

A = Ac.*fc + A_surf.*(1-fc);
% For when have input of lots of columns :-
%A(isnan(A))=A_surf;  %set NaN values to be clear-sky i.e. assuming no cloud)
inds_thresh = find(A<A_surf);
if length(A_surf)==1
    inds_A_surf=1;
else
   inds_A_surf = inds_thresh;
end
    
A(inds_thresh) = A_surf(inds_A_surf); %also don't allow below A_surf (the surface albedo)
        %(for cases when are providing columns with CF=1 for each one).

SW_out = SW_in .* A .*trans;



