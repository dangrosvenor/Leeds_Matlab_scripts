function [Ac_out,tau_out,A_out,SW_out,SW_out_clear_sky] = calc_SW2(fc,Wc,Nc,SW_in,A_surf,trans,cf_min,iprovide_tau,tau,i_Liu)
% Calculates the SW upwelling at TOA based (SW_out) on the cloud properties
% and SW downwelling at TOA (SW_in) and atmospheric transmissivity.
% NOW INCLUDING above cloud scattering due to 1 minus atmos transmissivity 
% fc is the cloud fraction in the range of 0 to 1.
% Nd in cm3
% W in g/m2
% SW_in - downward SW at TOA in W/m2
% N.B. making assumptions about cloud temperature and pressure in
% albedo_cloudy_func_Seinfeld function
% trans is the estimated transmission ratio of the atmosphere (i.e.
% SW_down_surf / SW_down_TOA )
% A_surf is the surface albedo.


%Call multiple times to account for multiple scattering between the surface
%and atmosphere
N=0;
SW_out = 0;
SW_out_clear_sky = 0;
SW_in_new = SW_in;
SW_single = 1e9;
SW_up_TOA_clear_sky_single = 1e9;
while maxALL(SW_single) > 0.2 | maxALL(SW_up_TOA_clear_sky_single) > 0.2 %for i=1:N
    N=N+1;
    [SW_in_new,Ac_out,tau_out,A_out,SW_single,SW_up_TOA_clear_sky_single] = SW_calc2_single(fc,Wc,Nc,SW_in_new,A_surf,trans,cf_min,iprovide_tau,tau,i_Liu);
    SW_out = SW_out + SW_single;   
    SW_out_clear_sky = SW_out_clear_sky + SW_up_TOA_clear_sky_single; %adding these in the same way as the all-sky fluxes    
end

A_out = SW_out ./ SW_in;
'';

function [SW_down_TOA_from_cloud_or_surf,Ac,tau,A,SW_up_TOA_total,SW_up_TOA_clear_sky] = SW_calc2_single(fc,Wc,Nc,SW_in,A_surf,trans,cf_min,iprovide_tau,tau,i_Liu)



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

%A = (Ac + (1-Ac).^2.*A_surf).*fc + A_surf.*(1-fc);
%A = (Ac ).*fc + A_surf.*(1-fc);
%The (1-Ac).^2.*A_surf term is to account for light that passes through the
%cloud, hits the surface and goes back through the cloud. So, now don't
%need the stuff below where don't allow the overall albedo to be less than
%the surface albedo.



% For when have input of lots of columns :-
%A(isnan(A))=A_surf;  %set NaN values to be clear-sky i.e. assuming no cloud)
% inds_thresh = find(A<A_surf);
% if length(A_surf)==1
%     inds_A_surf=1;
% else
%    inds_A_surf = inds_thresh;
% end
%     
% A(inds_thresh) = A_surf(inds_A_surf); %also don't allow below A_surf (the surface albedo)
%         %(for cases when are providing columns with CF=1 for each one).

        
% SW_surf_up_from_cloud_or_surf = SW_in .*trans .* A; %Outgoing. *trans since SW_in is TOA
% 
% 
% SW_up_single = SW_scatter + SW_surf_up_from_cloud_or_surf.*trans; %multiply by trans again for atmos losses
% 
% SW_down_TOA_from_cloud_or_surf = SW_surf_up_from_cloud_or_surf.*(1-trans); %Light emanting from the surface and scattered back down

%Radiation scattered directly back to TOA before reaching surface or cloud
%top - had neglected this before...
% This shouldn't matter if it's clear-sky or cloudy below
fabs = 0.3; %absorbed fraction
SW_scatter = SW_in .* (1 - trans) .* (1 - fabs);
SW_abs = SW_in .* (1 - trans) .* (fabs);

%Upwards flux above surface in clear-sky. 
SW_up_above_surf_clear_sky = SW_in .*trans .* A_surf .* (1-fc);

%Also account for light passing through clodu to hit the surface -
%assuming no clear-sky losses here.



N=0;
SW_down_below_cloud = SW_in .*trans .* (1-Ac) .* fc;
SW_up_above_cloud_multiple = SW_in .*trans .* (Ac) .* fc; %Outgoing. *trans since SW_in is TOA
while maxALL(SW_down_below_cloud)>0.2 %for i=1:N
    SW_up_above_cloud_multiple = SW_up_above_cloud_multiple + SW_down_below_cloud.*A_surf.*(1-fabs).*(1-Ac);
    SW_down_below_cloud = SW_down_below_cloud.*A_surf.*(1-fabs).*(Ac); 
    N=N+1;
end

SW_up_TOA_total = SW_scatter + (SW_up_above_surf_clear_sky + SW_up_above_cloud_multiple).*trans; %multiply by trans again for atmos losses
SW_down_TOA_from_cloud_or_surf = (SW_up_above_surf_clear_sky + SW_up_above_cloud_multiple) .*(1-trans); %Light emanating from the surface and scattered back down

%The flux if there was no cloud for comparison to model clear-sky fluxes.
SW_up_TOA_clear_sky = SW_scatter + (SW_up_above_surf_clear_sky ./ (1-fc) ).*trans .* (1 - fabs); %multiply by trans again for atmos losses
%Divide SW_up_above_surf_clear_sky by (1-fc) since is multiplied by this
%above and if fc=0 then 1-fc=1, so need to divide by (1-fc) and multiply by
%1 to make it equivalent to treated fc as 1.

A = SW_up_TOA_total ./ SW_in;
'';

