function [P_rad,qv_rad,T_rad,LWC_rad,LWC_rad_kg] = reverse_moist_adiabatic_path(P0,T0,Q0,LWC)
%function [P_rad,qv_rad,T_rad,LWC_rad,LWC_rad_kg] = reverse_moist_adiabatic_path(P0,T0,Q0,LWC)
%calculates a reverse moist adiabatic until all the LWC has evaporated. Then follows
%a dry adiabat. The starting conditions are specified by T0(K), P0(Pa), Q0(vapour MR in kg/kg)
%and starting LWC (g/m^3)
        
        dp=0.1e2;
        P_rad=[P0-dp:dp:1100e2]; %follow path to below sea level
        
        [LWC_rad,T_rad]=adLWC_PaulLawson_simple(P0,T0,P_rad); %solve the moist adiabat for the descent
        %note, produces negative LWC, which indicates evaporation of the LWC
        %this function gives temperature in degrees C
        
        %%% need to convert LWC into kg/kg
            %first for the initial point
            rho_eq=density(P0,T0);
            LWCeq_kg=LWC/1000 /rho_eq;

            %now for the points in the descent
            rho_ad=density(P_rad,T_rad+273.15);
            LWC_rad_kg = LWC_rad/1000 ./ rho_ad;
        
        %vapour MR during the descent = starting vapour minus the LWC produced (LWC negative in this case)    
        qv_rad = Q0-LWC_rad_kg; 
        
        %look for when the combination of the initial LWC and the LWC lower down
        %is zero (i.e. all LWC evaporated)
        [minLWC,ii2]=min(abs(LWC_rad_kg+LWCeq_kg(1))); 
        if abs(minLWC)>0.1e-3
            fprintf(1,'\n**** WARNING - not all LWC evaporated on moist descent. %f left ***\n',minLWC);
            fprintf(1,'**** aborting. Override if required, or extend to a higher pressure ***\n');
            return
        end
        
        %so from index ii2 onwards we have are following the dry adiabat
        pot0 = (T_rad(ii2)+273.15) * (1000e2/P_rad(ii2))^0.286;
        T_rad(ii2+1:end) = pot0./ (1000e2./P_rad(ii2+1:end)).^0.286 - 273.15;
%        rho_cb = density(P_rad(ii2),T_rad(ii2));
        qv_rad(ii2+1:end) = qv_rad(ii2);
        
        
        