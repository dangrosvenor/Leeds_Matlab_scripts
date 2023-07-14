clear clims_delta

switch var_ukesm
        case 'Nd_cf_weighted_UKESM_ztop'
            %load('/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm.mat');
        case {'Nd_cf_weighted_UKESM','scldncl','Nd_clw_weighted_ESGF','Nd_clw_weighted_ESGF_no_dz','Nd_clw_weighted_ESGF_no_dz_no_ice_total_column_to_zdomain_top'}
            var_str = 'N_d';
            units_str_trend = 'cm^{-3} yr^{-1}';
            units_str = 'cm^{-3}';
            switch yr_start_trend
                case 2003
                    clims=[-6 6];
                case {1850,1851}
                    switch yr_end_trend
                        case 2014
                            clims=[-1.5 1.5];
                        case 1970
                            clims=[-2 2];
                    end
            end
            
        case {'reffclwtop'}
            var_str = 'Liquid cloud top r_e';
            units_str_trend = '\mum yr^{-1}';
            units_str = '\mum';
            clims=[-6 6];               

        case 'calipso_low_cloud_amount'
            %var_str = 'COSP CALIPSO low cloud fraction';
            var_str = 'Low cloud fraction';
            units_str_trend = 'yr^{-1}';
            units_str = '';
            clims=[-0.003 0.003];
            
     case {'calipso_total_cloud_amount','clt'}
         %var_str = 'COSP CALIPSO low cloud fraction';
         var_str = 'Total cloud fraction';
         var_str = 'f_{c}';
         units_str_trend = 'yr^{-1}';
         units_str = '';
         clims=[-0.003 0.003];
            
     case {'SW_up_TOA','rsut','SWTOA Calc'}
         %var_str = 'F_{SW TOA upwelling}';
         var_str = 'F_{SW \uparrow}';
         units_str_trend = 'W m^{-2} yr^{-1}';
         units_str = 'W m^{-2}';
         clims=[-0.7 0.7];
         clims=[-0.4 0.4]; 
         clims_delta = [-25 25];
         
     case {'rsutcs'}
         %var_str = 'Clear-sky SW TOA up';
         var_str = 'F^{clear-sky}_{SW\uparrow}';
         units_str_trend = 'W m^{-2} yr^{-1}';
         units_str = 'W m^{-2}';
         clims=[-0.7 0.7];
         clims=[-0.4 0.4];
         clims=[-0.2 0.2];
         clims=[-0.05 0.05];
         
      case {'rsds'}
         var_str = 'F_{SW\downarrow surf}';
         units_str_trend = 'W m^{-2} yr^{-1}';
         units_str = 'W m^{-2}';
         clims=[-0.7 0.7];
         clims=[-0.4 0.4];      
         
     case 'SO2_low_anthropogenic_emissions'
         var_str = 'SO2 emissions';
         units_str_trend = 'tonnes km^{-2} yr^{-1} yr^{-1}';
         units_str = 'tonnes km^{-2} yr^{-1}';
         clims = 1e-12*[-1 1]*3600*24*365*1e6/1e3; %convert from kg/m2/s to tonnes/km2/yr
         
         
     case 'ts'
         %var_str = 'Surface Temperature';
         var_str = 'T';
         units_str_trend = 'K yr^{-1}';
         units_str = 'K';
         clims = [-0.1 0.1];
         
         
         
     case 'dust_od550'
         var_str = 'Dust 550nm AOD';
         units_str_trend = 'yr^{-1}';
         units_str = '';
         clims = [-0.008 0.02];
         
     case 'lwp'
         %var_str = 'LWP';
         var_str = 'L_{all-sky}';
         units_str_trend = 'g m^{-2} yr^{-1}';
         units_str = 'g m^{-2}';
         clims = [-1 1];
         
     case 'lwpic'
         %var_str = 'LWPic';
         var_str = 'L';
         units_str_trend = 'g m^{-2} yr^{-1}';
         units_str = 'g m^{-2}';
         clims = [-1 1];    
         
     case 'clwvi'
         var_str = 'LWP_{clwvi}';
         units_str_trend = 'g m^{-2} yr^{-1}';
         units_str = 'g m^{-2}';
         clims = [-1 1];
         
     case 'prw'
         var_str = 'Water Vapour Path';
         units_str_trend = 'kg m^{-2} yr^{-1}';
         units_str = 'kg m^{-2}';
         clims = [-1 1];        
         
     case 'od550aer'
         var_str = 'Aerosol Optical Depth at 550nm';
         units_str_trend = 'yr^{-1}';
         units_str = '';
         clims = [-1 1];      
         
    case 'od550tot'
        %var_str = 'Aerosol+Dust Optical Depth,550nm';
        var_str = '\tau_{a}';
        units_str_trend = 'yr^{-1}';
        units_str = '';
        clims = [-1 1];
         
     case 'od550aerso'
         var_str = 'Stratospheric AOD at 550nm';
         units_str_trend = 'yr^{-1}';
         units_str = '';
         clims = [-1 1];
         
     otherwise
         var_str = 'NOT SET';
         units_str_trend = '';
         units_str = '';
         clims = [-1 1];
           
end
    


if ~exist('clims_delta')
   clims_delta = clims*30; 
end