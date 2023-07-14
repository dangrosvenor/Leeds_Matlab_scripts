function [Ac,tau,A_out,SW_out,SW_out_clear_sky] = calc_SW3(opts,aod,aaod,f_upscatter,trans_clear_sky,fc,Wc,Nc,SW_in,A_surf,trans)
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
% A_surf is the surface albedo - usually only used for the multiple scatter
% below clouds (otherwise for clear-skies it uses the provided clear-sky
% SW).

%   %Convert all of the variable names in the input structure to actual names
%     %for ease of use
%     name_struc='optional'; %The name of the structure
%     names = eval(['fieldnames(' name_struc ');']);
%     for i=1:length(names)
%         eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
%         eval(eval_str);
%     end


if ~isfield(opts,'cf_min')
    opts.cf_min = 0.01;
end
if ~isfield(opts,'iprovide_tau')
    opts.iprovide_tau = 0;
end
if ~isfield(opts,'i_Liu')
    opts.i_Liu=0;
end
if ~isfield(opts,'iprovide_SWclear')
    opts.iprovide_SWclear = 0;
end
if ~isfield(opts,'idiv_cosSZA')
    opts.idiv_cosSZA = 0;
end

Nc = Nc*1e6; %convert to per m3
Wc = Wc/1e3; %convert to kg/m2

if opts.iprovide_tau==0
    [Ac,tau] = albedo_cloudy_func_Seinfeld(Wc,Nc,opts.i_Liu);
else
    Ac = opts.tau ./(opts.tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)
end

tau(fc<opts.cf_min)=NaN;
Ac(fc<opts.cf_min)=0; %When there is no cloud the cloud albedo will be NaN, so set cloud albedo to zero to return the clear-sky value


iexperimental_multiple_scatter=0; %Seems to work quite well actually
if iexperimental_multiple_scatter==1
   SSA_cloud = 1;
   f_up_cloud = 1;
   r = Ac.*SSA_cloud.*f_up_cloud; %Global mean values of upscatter start at 0.5 for small particles and then decrease. 0.4 might be reasonable.
    %Transmissivity of aerosol layer (single layer) :-
    t = (1-Ac) + SSA_cloud.*(1-f_up_cloud).* Ac;
    %Mulitple scattering between aerosol layer and surface, plus
    %atmospheric transmissivity effect (i.e., even with no aerosol)
    Ac =  r + (t.^2.*A_surf)./(1-A_surf.*r); 
end

%% Clear-sky aerosol induced scatter following Seinfeld and Pandis textbook Chapter 24

if opts.iprovide_SWclear==1    
    SW_out_clear_sky = opts.SWclear; %But can just get this from the clear sky surface flux of the model
else
    
    %Eqn 24.4
    SSA = (aod - aaod ) ./ aod; %from Eqn 1 of Lacagnina JGR 2015; doi:10.1002/2015JD023501
    SSA(aod<0.001)=1;
    eAOD = exp(-aod);
    %Reflectivity of aerosol layer (single pass) :-
    r = (1 - eAOD).*SSA.*f_upscatter; %Global mean values of upscatter start at 0.5 for small particles and then decrease. 0.4 might be reasonable.
    %Transmissivity of aerosol layer (single layer) :-
    t = eAOD + SSA.*(1-f_upscatter).*(1 - eAOD);
    %Mulitple scattering between aerosol layer and surface, plus
    %atmospheric transmissivity effect (i.e., even with no aerosol)
    SW_out_clear_sky = trans_clear_sky.^2 .* ( r + (t.^2.*A_surf)./(1-A_surf.*r)).*SW_in;
    if opts.idiv_cosSZA == 1
        SW_out_clear_sky = SW_out_clear_sky ./ opts.cosSZA;
    end
end

%% Combined upwards flux

SW_out = Ac .* fc .* SW_in .* trans.^2 + (1-fc).*SW_out_clear_sky;
A_out = SW_out ./ SW_in;




