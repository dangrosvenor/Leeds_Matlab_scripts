

if isstr(time_round)
    time_round2=time_round;
else
    time_round2=datestr(time_round);
end
    


%%
switch var_UM
    case 'wind speed';    

        U_glm = nc{'U'}(it,:,:);
        V_glm = ncV{'V'}(it,:,:);
        V_glm = ncV{'V'}(it,1:end-1,:);

        dat_modis = sqrt(U_glm.^2 + V_glm.^2);
        
        titlenam_driver = ['Wind speed (m s^{-1}) at ' time_round2 time_format_str];
        units_str_plot = '';

    case 'accum_number_ukca'
        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['Soluble accumulation mode number (# per #molecules air) at ' time_round2 time_format_str];
        
        clim_min = 0;
        clim_max = 3e-17;
        
    case 'accum_number_ukca diff'
        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['Difference']
        
        clim_min = -0.5e-17;
        clim_max = 0.5e-17;        
        
    case 'accum_mass_H2SO4_ukca_xxx'
        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['Soluble accumulation mode H2SO4 MMR (kg kg^{-1}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 5e-9;
        
    case 'accum_mass_OC_ukca'
        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['Soluble accumulation mode OC MMR (kg kg^{-1}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 2e-9; 
        
    case 'accum_mass_BC_ukca'
        
        %dat_modis = nc{var_UM}(it,:,:);
        clim_min = 0;
        clim_max = 5e-10;

        titlenam_driver = ['Soluble accumulation mode BC MMR (kg kg^{-1}) at ' time_round2 time_format_str];

    case 'emissions_BC_biomass'    
        dat_modis = BC_model_level01_time08;
        titlenam_driver = ['BC biomass emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];
        clim_min = 0;
        clim_max = 5e-11;
        
        
    case 'emissions_BC_biofuel'
        
        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 5e-12; 
        
    case 'emissions_BC_fossil'

        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['BC fossil fuel emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 5e-12;

    case 'emissions_DMS'

        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['DMS emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 5e-12;

    case 'emissions_Monoterp'

        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['Monoterpene emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 5e-12;

    case 'emissions_OC_biofuel'

        %dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 5e-12;

%     case 'emissions_BC_biofuel'
% 
%         %dat_modis = nc{var_UM}(it,:,:);
%         titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];
% 
%         clim_min = 0;
%         clim_max = 5e-12;
% 
%     case 'emissions_BC_biofuel'
% 
%         %dat_modis = nc{var_UM}(it,:,:);
%         titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' time_round2 time_format_str];
% 
%         clim_min = 0;
%         clim_max = 5e-12;
        
        
    case 'LWP'

        %%dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['LWP (g m^{-2}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 300;
        
    case 'LWP diff'

        %%dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['\DeltaLWP (g m^{-2}) at ' time_round2 time_format_str];

        clim_min = -100;
        clim_max = 100;  
        
        
    case 'SW_down_surf'

        %%dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['SW down at surface (W m^{-2}) at ' time_round2 time_format_str];

        clim_min = -30;
        clim_max = 10;  
        
   case 'SW_down_surf diff'

        %%dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['\Delta SW down surf (W m^{-2})'];

        clim_min = -400;
        clim_max = 400;  


    case 'Nd_lwc_weighted_UKCA'

        %%dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['N_d (cm^{-3}) at ' time_round2 time_format_str];

        clim_min = 0;
        clim_max = 300;
        
    case 'Nd_lwc_weighted_UKCA diff'

        %%dat_modis = nc{var_UM}(it,:,:);
        titlenam_driver = ['\Delta N_d (cm^{-3}) (model minus MODIS) at ' time_round2 time_format_str];

        clim_min = -100;
        clim_max = 100;
        
    case 'SO2_column'
        dat_modis = dat_modis *1000 / 0.0286; %convert to Dobson Units 
        % 	1 DU = 2.69e20 molecules m^-2, molecular mass of SO2 = 64.066 g/mol.
        %   So, 1 DU for SO2 is equivalent to a column mass of :-m = 2.687e20/6.023e23 *64.066 = 0.0286 g/m2)        
        titlenam_driver = ['SO_2 column integrated mass (DU) at ' time_round2 time_format_str];
        
        clim_min = 0;
        clim_max = 2.5e-4*1000 / 0.0286;
        
    case 'accum_mass'
        titlenam_driver = ['CASIM accumulation mode column integrated mass (kg m^{-2}) at ' time_round2 time_format_str];
        clim_min = 0;
        clim_max = 4e4;
        
    case 'accum_mass_H2SO4_ukca'
        titlenam_driver = ['UKCA accumulation mode H_2SO_4 soluble column integrated mass (kg m^{-2}) at ' time_round2 time_format_str];
        clim_min = 0;
        clim_max = 4e-5;   
        
    case 'LS_surf_rain_rate'

        dat_modis = dat_modis * 3600; %convert from mm per s to mm per hour
        titlenam_driver = ['LS surface rain rate (kg m^{-2} hr^{-1}) at ' time_round2 time_format_str]; %same as mm/hr

        clim_min = 0;
        clim_max = 300;
      
    otherwise
        titlenam_driver = remove_character([var_UM ' at ' time_round2 time_format_str],'_',' ');
        
end




