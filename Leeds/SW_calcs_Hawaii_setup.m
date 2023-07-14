
%% Atmospheric transmissivity
switch tr_for_SW_bias
    case 'constant value'
        tr_atmos_mean = tr_atmos_constant_val; %if constant value
    case 'time mean from model'
        SW_thresh=200;
        %transmission_atmos = SW_clean_clear_PD_ALL ./ SW_in_PD_ALL;
        %transmission_atmos(SW_in_PD_ALL<SW_thresh)=0.5;
        %tr_atmos_mean = meanNoNan(transmission_atmos,3);
        use model clear-sky
    case {'time varying monthly','time varying annual mean'}
        transmission_atmos_up = sw_TOA_up_clear_dat_ens_mean ./ sw_surf_up_clear_dat_ens_mean;
        %Upwelling from at TOA divided by upwelling at surface in
        %clear-skies.
        
        transmission_atmos_down = sw_surf_down_clear_dat_ens_mean ./ SW_down_TOA_dat_ens_mean;
        %Downwelling to the surface divided by downwelling at TOA in
        %clear-skies.
        
        %So, will get two measures of transmissivity here. The TOA
        %calculation uses a cloud top SW downwelling, so could just use
        %the surface downwelling for this (although may underestimate a
        %little).
        %Then can provide the upwelling trans value to the SW calc
        
        %Actually, can only really use the downwards transmission since upwards includes
        %the scattered light from the sun and not just that originating
        %from the surface.
end

%% Surface albedo
switch albedo_type
    case 'constant'
        % Also need the clear-sky albedo - can estimate from cloud free regions?
        surf_albedo = 0.15; %Guess for now - value for CF quite sensitive to this (ranging up to 30%), so should prob
        %check what this is from the model
        surf_albedo = 0.12; %value estimated from the 12 noon snapshot from VOCALS reginal run - N.B. - this uses the estimated transmission since
        %it is the surface albedo assuming no atmosphere above
        %surf_albedo = 0.0;
        surf_albedo = 0.09;
        
    case 'calculated all time steps'
        %Estimate of surface albedo using surface fluxes in clear-sky
        surf_albedo2  = SW_up_clean_clear_surf_PI_ALL ./ SW_down_clean_clear_surf_PI_ALL;
        inan = find(SW_down_clean_clear_surf_PI_ALL<200);
        surf_albedo2(inan) = NaN;
        %Don't want NaNs - so do a time mean and replicate for all times
        me = meanNoNan(surf_albedo2,3);
        surf_albedo = repmat(me,[1 1 size(surf_albedo2,3)]);
        
        
    case 'calculated annual mean'
        surf_albedo_monthly  = sw_surf_up_clear_dat_ens_mean ./ sw_surf_down_clear_dat_ens_mean;
        
        clear surf_albedo_ann
        for iy=1:length(years_sw_calc)
            istart=(iy-1)*12+1;
            %surf_albedo_ann(iy,:,:) = meanNoNan(surf_albedo_monthly(istart:istart+11,:,:),1,'mean',1,1,SW_down_TOA_dat_ens_mean(istart:istart+11,:,:));
            surf_albedo_ann(iy,:,:) = min(surf_albedo_monthly(istart:istart+11,:,:),[],1);
        end
        %trans_orig = repmat(trans_orig_ann,[1 1 1 12]);
        surf_albedo = NaN*ones(size(surf_albedo_monthly));
        for im=1:size(surf_albedo_monthly,1)
            iy = ceil(im/12);
            surf_albedo(im,:,:) = surf_albedo_ann(iy,:,:);
        end
        
        
        
end


