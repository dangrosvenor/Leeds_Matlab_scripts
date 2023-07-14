function [P_rad,qv_rad,T_rad,LWC_rad] = forward_dry_then_moist_adiabatic_path(P0,T0,Q0,Pend)
%function [P_rad,qv_rad,T_rad,LWC_rad,LWC_rad_kg] = reverse_moist_adiabatic_path(P0,T0,Q0,LWC)
%calculates a reverse moist adiabatic until all the LWC has evaporated. Then follows
%a dry adiabat. The starting conditions are specified by T0(K), P0(Pa), Q0(vapour MR in kg/kg)
%and starting LWC (g/m^3)
f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

        dp=-0.1e2;
        P_rad=[P0:dp:Pend]; %follow path up to Pend

        
    %following the dry adiabat first
        pot0 = T0 * (1000e2/P0)^0.286;
        T_rad = pot0./ (1000e2./P_rad).^0.286;
%        rho_cb = density(P_rad(ii2),T_rad(ii2));
        qv_rad = Q0*ones([1 length(T_rad)]);

        qsat = SatVapPress(T_rad,'goff','liq',P_rad,1)/f;
        [min_val,imin]=min(abs(qv_rad-qsat));        
        
        LWC_rad = zeros([1 length(T_rad)]);
        
        T_rad=T_rad-273.15;
        [LWC_rad(imin:end),T_rad(imin:end)]=adLWC_PaulLawson_simple(P_rad(imin),T_rad(imin)+273.15,P_rad(imin:end)); %solve the moist adiabat for the descent
        %note, produces negative LWC, which indicates evaporation of the LWC
        %this function gives temperature in degrees C
        qv_rad(imin:end)=qsat(imin:end);  
        
        