function []=ACSIS_dat_trends_load_ensemble_esgf(var_UM,MIP,expt)
%Run from ACSIS_dat_trends_load_ESGF_ensemble_multi_vars.m

output_period = 'all';
%output_period = 'recent';
%output_period = 'SW_down_TOA_partial';

i_annual_data=0; %(setting default). Flag to say that the data is already annual averages (not monthly).

%Loading and scaling options
opts.isort_time=1;
iuse_mat_file=0;
fscale_fac = 1;
opts.lat_var = 'lat';
opts.lon_var = 'lon';
    
    
    
switch MIP
    case 'CMIP6'                
        
        ens_runs={'r1i1p1f2','r3i1p1f2','r5i1p1f3','r7i1p1f3','r9i1p1f2'...
            'r2i1p1f2','r4i1p1f2','r6i1p1f3','r8i1p1f2','r11i1p1f2','r16i1p1f2'...
            ,'r18i1p1f2','r10i1p1f2','r12i1p1f2','r17i1p1f2','r19i1p1f2'}; ens_str='';
        
        %ens_runs={'r1i1p1f2'}; ens_str='';
        
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/'];
        %dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MOHC/UKESM1-0-LL/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;
         
    case 'CMIP6 HadGEM3-GC31-LL'
        cmip6_str = 'HadGEM3-GC31-LL';
        
         ens_runs={'r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3','r5i1p1f3'}; ens_str='';
               
        
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MOHC/' cmip6_str '/historical/'];        
        %Directory to save to        
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];        
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;         
         
    case 'CMIP6 UKESM'
        cmip6_str = 'UKESM1-0-LL';
        
        ens_runs={'r1i1p1f2','r3i1p1f2','r5i1p1f3','r7i1p1f3','r9i1p1f2'...
            'r2i1p1f2','r4i1p1f2','r6i1p1f3','r8i1p1f2','r11i1p1f2','r16i1p1f2'...
            ,'r18i1p1f2','r10i1p1f2','r12i1p1f2','r17i1p1f2','r19i1p1f2'}; ens_str='';
        
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MOHC/' cmip6_str '/historical/'];        
        %Directory to save to        
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];        
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;  
         
    
        
        
        
        
case 'CMIP6 UKESM1-AerChemMIP_control'
        cmip6_str = 'UKESM1-AerChemMIP_control';
        
        %ens_runs={'r1i1p1f2','r3i1p1f2','r5i1p1f3','r7i1p1f3','r9i1p1f2'...
        %    'r2i1p1f2','r4i1p1f2','r6i1p1f3','r8i1p1f2','r11i1p1f2','r16i1p1f2'...
        %    ,'r18i1p1f2','r10i1p1f2','r12i1p1f2','r17i1p1f2','r19i1p1f2'}; ens_str='';
        
        ens_runs={'r1i1p1f2','r2i1p1f2','r3i1p1f2'}; ens_str='';        
  
        %Directory to process :-
        switch var_UM
            case {'Nd_clw_weighted_ESGF','Nd_clw_weighted_ESGF_total_column_to_zdomain_top'}
                dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MOHC/UKESM1-0-LL/historical/']; 
            otherwise
                dir_data = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/'];
                
        end    
        
        
        
        %Directory to save to        
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];        
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;  
         
    case 'CMIP6 BCC-ESM1'
        cmip6_str = 'BCC-ESM1';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/BCC/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;     
         
    case 'CMIP6 CNRM-CM6-1'
        cmip6_str = 'CNRM-CM6-1';
        
        ens_runs={'r1i1p1f2','r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2','r6i1p1f2'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/CNRM-CERFACS/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;      
         
     case 'CMIP6 CNRM-ESM2-1'
        cmip6_str = 'CNRM-ESM2-1';
        
        ens_runs={'r1i1p1f2','r2i1p1f2','r3i1p1f2'}; ens_str='';
        
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/CNRM-CERFACS/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;    
         
   case 'CMIP6 EC-Earth3-AerChem'
        cmip6_str = 'EC-Earth3-AerChem';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/EC-Earth-Consortium/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;    
         
         
 case 'CMIP6 MPI-ESM-1-2-HAM'
        cmip6_str = 'MPI-ESM-1-2-HAM';
        
        ens_runs={'r1i1p1f1','r2i1p1f1','r3i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/HAMMOZ-Consortium/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;  
         
    case 'CMIP6 IPSL'
        cmip6_str = 'IPSL-CM6A-LR';
        
        ens_runs={'r1i1p1f1'}; ens_str=''; %Are other members, but not processed for some reason
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/IPSL/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;     
         
    case 'CMIP6 MIROC6'
        cmip6_str = 'MIROC6';
        
        ens_runs={'r1i1p1f1','r10i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1'}; ens_str=''; 
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MIROC/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;  
  
    case 'CMIP6 MIROC-ES2L'
        cmip6_str = 'MIROC-ES2L';
        
        ens_runs={'r1i1p1f2'}; ens_str=''; %Are other members, but not processed for some reason
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MIROC/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014; 
         
    
    case 'CMIP6 MRI-ESM2-0';
        cmip6_str = 'MRI-ESM2-0';
        
        ens_runs={'r1i1p1f1'}; ens_str=''; %Are other members, but not processed for some reason
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MRI/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014; 
             
         
    case 'CMIP6 CESM2';
        cmip6_str = 'CESM2';
        
        ens_runs={'r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'}; ens_str=''; ens_str=''; %Are other members, but not processed for some reason
        %ens_runs={'r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'}; ens_str=''; ens_str=''; %Are other members, but not processed for some reason
        %ens_runs={'r1i1p1f1'}; ens_str=''; ens_str=''; %
        %ens_runs={'r1i1p1f1','r2i1p1f1','r3i1p1f1'};
       
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NCAR/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;   
         
    case 'CMIP6 CESM2-FV2';
        cmip6_str = 'CESM2-FV2';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NCAR/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
        start_year = 1850;
        end_year = 2014;
        
    case 'CMIP6 CESM2-WACCM';
        cmip6_str = 'CESM2-WACCM';
        
        %ens_runs={'r1i1p1f1','r2i1p1f1','r3i1p1f1'}; ens_str='';
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NCAR/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014; 
         
    case 'CMIP6 CESM2-WACCM-FV2';
        cmip6_str = 'CESM2-WACCM-FV2';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NCAR/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;          

    case 'CMIP6 NorESM2-MM';
        cmip6_str = 'NorESM2-MM';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NCC/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;          
         
    case 'CMIP6 NIMS-KMA UKESM1-0-LL';
        cmip6_str = 'UKESM1-0-LL';
        
        ens_runs={'r13i1p1f2','r14i1p1f2'}; ens_str='';
        ens_runs={'r13i1p1f2','r14i1p1f2','r15i1p1f2'}; ens_str='';
        ens_runs={'r14i1p1f2','r15i1p1f2'}; ens_str='';
        %ens_runs={'r13i1p1f2'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NIMS-KMA/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_NIMS-KMA' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;
         
    case 'CMIP6 GFDL-CM4';
        cmip6_str = 'GFDL-CM4';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NOAA-GFDL/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;   
         
    case 'CMIP6 GFDL-ESM4';
        cmip6_str = 'GFDL-ESM4';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/NOAA-GFDL/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014; 
         
    case 'CMIP6 INM-CM4-8'
        cmip6_str = 'INM-CM4-8';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/INM/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014; 

    case 'CMIP6 INM-CM5-0'
        
        cmip6_str = 'INM-CM5-0';
        
        ens_runs={'r1i1p1f1'}; ens_str='';
        
        %Directory to process :-        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/INM/' cmip6_str '/historical/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' cmip6_str '_' output_period ens_str];                        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014; 

         
    case 'DAMIP'               
        
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/UM/UKESM/DAMIP/' expt '/'];  
        
        ens_runs={'r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'}; ens_str='';
        
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' MIP '_' expt '_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
        start_year = 1850;
        end_year = 2014;  
        
        
    case 'AerChemMIP'               
        
        %Directory to process :-
        
        dir_data = ['/home/disk/eos15/d.grosvenor/ESGF/MOHC/UKESM1-0-LL/' expt '/'];          
        
        ens_runs={'r1i1p1f2','r2i1p1f2','r3i1p1f2'}; ens_str='';
        
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' MIP '_' expt '_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
        start_year = 1850;
        end_year = 2014;          
         
         
    case 'HADGEM3_GC31_LL'
                
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/UM/HadGEM-GC31/HadGEM-GC31-LL/historical/'];  
        
        %ens_runs={'r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3','r5i1p1f3'}; ens_str='';
        ens_runs={'r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'}; ens_str='';
        
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' MIP '_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
         start_year = 1850;
         end_year = 2014;           
        
    case 'PI'             
        %PI control - just need to calculate the time mean and std. dev (only one
        %member). Runs from 1960 to 2709 (actaully year 3839 for ESGF it seems) giving 750 years, or 9000 months.
        
        ens_runs={'r1i1p1f2'}; ens_str='_PI_control'; 
        
        %ESGF data :-
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_PIcontrol/'];
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_' ens_str];
        
        load_type = 'individual netCDF files'; opts.cat_dim=0; %set to zero since each file contains just latxlon data.
            %This sets the time dim as the first one and cats all files
            %(one time per file here - if had more than one time per file
            %then would set cat_dim to the time dimension.
                
        i_annual_data=1; %Flag to say that the data is already annual averages (not monthly).
        
       
        % Data from .pp files :-
        %ens_runs={'picontrol_r1i1p1f2_u-aw310'}; ens_str='_PI_control'; output_period = 'all_PI';
        %PI data located in /home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/UKESM/picontrol_r1i1p1f2_u-aw310
        
        start_year = 1960;
        end_year = 3839;
        
        
        
        
    case 'AMIP'
        dir_data = ['/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_AMIP_data/'];
        
        ens_runs={'1','dummy'};
        
        output_period = 'all';
        %output_period = 'recent';
        %output_period = 'SW_down_TOA_partial';
        
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_AMIP_' output_period];
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        
        start_year = 1979;
        end_year = 2014;
        
    case 'Nudged'               
        
        %Directory to process :-
        dir_data = ['/home/disk/eos15/d.grosvenor/UM/UKESM/UKESM_ACSIS_nudged_run_James_Keeble/' expt '/'];        
        
        
        ens_runs={'1','dummy'}; ens_str='';
        
        %Directory to save to
        savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' MIP '_' expt '_' output_period ens_str];
                
        
        
        load_type = 'mat';
        load_type = 'merged netCDF';
        %load_type = 'individual netCDF files'; opts.cat_dim=1;
        
        start_year = 1982; 
        end_year = 2014;   %full run
        %end_year = 2001; %currently only got this far for u-by844
        
        opts.lat_var = 'Latitude';
        opts.lon_var = 'Longitude';
        
end

ind_start_offset = 0; %offset for first time index of data (since SW TOA starts with a partial year)

%Whether is a single ensemble member process job
if length(ens_runs)>1
    i_single_ens=0;
else
    i_single_ens=1;
end

savefile = [savefile_pre_str '_' var_UM '.mat'];
for iens=1:length(ens_runs)
    
    pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
    switch MIP
        case 'CMIP6'
            dirUM = [dir_data ens_runs{iens}];
           
        case {'AMIP','Nudged'}
            dirUM = [dir_data];
            
        case 'PI'
            dirUM = [dir_data ens_runs{iens}];
            
        case 'DAMIP'
            dirUM = [dir_data ens_runs{iens}];
            
        otherwise
            dirUM = [dir_data ens_runs{iens} '/output/'];
            dirUM = [dir_data ens_runs{iens}];
            %error('Need to set something here!')
                   
    end
    
    
    
    

    var_UM_load = var_UM; %Might get reset below.
    switch var_UM
        case {'prra','tos'}
            opts.lat_var = 'latitude';
            opts.lon_var = 'longitude';
            opts.fill_range = [1e19 1e20];
            opts.grid_data_to_UKESM=1;
            
        case {'reffclwtop','cldnvi'}
            opts.fill_range = [1e19 1e21];
            
        case 'scldncl'
            iuse_mat_file=1;
            mat_file = [dirUM '/scldncl/monthly_Nd.mat'];
            fscale_fac = 1e-6; %convert from per m3 to per cm3
            
        case {'lwpic'}
            fscale_fac = 1e3; %convert from kg/m2 to per g/m2
            var_UM_load = 'lwp';
            
        case {'clwvi','lwp'}
            fscale_fac = 1e3; %convert from kg/m2 to per g/m2
        case {'Nd_clw_weighted_ESGF','Nd_clw_weighted_ESGF_total_column_to_zdomain_top','Nd_clw_weighted_ESGF_no_dz','Nd_clw_weighted_ESGF_no_dz_div_CF','Nd_clw_weighted_ESGF_no_dz_div_CF_total_column_to_zdomain_top',...
                'Nd_clw_weighted_ESGF_no_dz_div_CF_no_ice_total_column_to_zdomain_top','Nd_clw_weighted_ESGF_no_dz_no_ice_total_column_to_zdomain_top'}
            %.nc files created on Jasmin using Python
            opts.lat_var = 'Latitude';
            opts.lon_var = 'Longitude';
            
        case {'od550tot'}
            var_UM_load = 'od550aer';
            
    end
    
    
    if iuse_mat_file==1
        load(mat_file,'dat_global');
    else
        opts.time_var='time';
        opts.time_ref = datenum('01-Jan-1850');
        opts.time_fconv = 1; %conversion multiplier to get to days
        dat_global = UM_load_merged_netCDF(dirUM,var_UM_load,run_type,load_type,[],[],var_UM_load,opts); %data is ordered [time lat lon]. 180 times (monthly over 15 years)
    end
    
    if iens==1
        dat_ens = NaN*ones([length(ens_runs) size(dat_global.dat)]);
    end
    
    %deal with the ones that need masks (cloud fraction ,etc)
    switch var_UM
        case {'calipso_low_cloud_amount','calipso_mid_cloud_amount','calipso_high_cloud_amount','calipso_total_cloud_amount'}
            var_UM_mask = [var_UM '_mask'];
            dat_mask = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type,[],[],var_UM_mask,opts);
            %store all the dat data in a big array
            %dat_ens(iens,:,:,:) = fscale_fac .* dat_global.dat ./ dat_mask.dat;
            dat_ens(iens,:) = fscale_fac .* dat_global.dat(:) ./ dat_mask.dat(:);
        case {'lwpic','clwvi_ic'}
            var_UM_cf = ['clt'];
            dat_cf = UM_load_merged_netCDF(dirUM,var_UM_cf,run_type,load_type,[],[],var_UM_cf,opts);
            %total CF in percent.
            %store all the dat data in a big array
            imin = find(dat_cf.dat < 2); %limit min CF
            dat_cf.dat(imin) = 2;
            %dat_ens(iens,:,:,:) = fscale_fac .* dat_global.dat ./ (dat_cf.dat * 0.01);
            dat_ens(iens,:) = fscale_fac .* dat_global.dat(:) ./ (dat_cf.dat(:) * 0.01);
            %Divide by the total cloud fraction to get the in-cloud LWP
            
        case {'od550tot'} %add the normal aod and the dust aod
            var_UM2 = ['od550dust'];
            dat2 = UM_load_merged_netCDF(dirUM,var_UM2,run_type,load_type,[],[],var_UM2,opts);
            dat_ens(iens,:) = fscale_fac .* (dat_global.dat(:) + dat2.dat(:));
            
        otherwise
            %store all the dat data in a big array
            %dat_ens(iens,:,:,:) = fscale_fac .* dat_global.dat;
            dat_ens(iens,:) = fscale_fac .* dat_global.dat(:);
    end
    
end

if i_single_ens==0
    %Calculate ensemble mean and the inter-ensemble std dev
    [dat_ens_mean,dat_ens_Ndatap dat_ens_std] = meanNoNan(dat_ens,1);
else
    dat_ens_mean = squeeze(dat_ens);
end


%% Calculate annual means, etc.
gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
%Had this the wrong way around :-
%[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
[gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);

%manually set the time since matlab time doesn't work for the 360 day
%calendar (I think this is why the times are weird).
years_ukesm = repmat([start_year:end_year],[12 1]);
years_ukesm = years_ukesm';
years_ukesm_1d = years_ukesm(:,1);
months_ukesm = repmat([1:12],[length(years_ukesm) 1]);
days_ukesm = ones(size(months_ukesm));
time_ukesm = datenum(years_ukesm,months_ukesm,days_ukesm);

%For ISCCP pressure vs tau cloud fraction data
%size(dat_ens) = [nens nmonths ntau plev nlat nlon]

sdat=size(dat_ens);

%dat_std_ukesm = NaN*ones([length(years_ukesm_1d) size(dat_ens_mean,2) size(dat_ens_mean,3)]); %not used

if i_annual_data==1
    dat_annual = dat_ens_mean;
    dat_annual_ens = dat_ens;
    dat_ukesm = dat_ens_mean;
else
    dat_annual = NaN*ones([length(years_ukesm_1d) sdat(3:end)]);
    dat_ukesm = NaN*ones([length(years_ukesm_1d) 12 sdat(3:end)]);
    dat_annual_ens = NaN*ones([sdat(1) length(years_ukesm_1d) sdat(3:end)]);
    
    for iy=1:length(years_ukesm_1d) %Loop over years and use indices for the monthly data to create annual averages
        %for im=1:12
        tind_01 = ind_start_offset + (iy-1)*12 + 1;
        tind_02 = tind_01+11;
        %dat_annual(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %%calculate the annual means keeping each ensemble separate
        %dat_annual_ens(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        %Replaced the above with just one : at the end to allow for cases where
        %have more than just lat and lon dimensions (e.g., pressure, optical
        %depth for ISCCP). Should still work the same.
        dat_annual(iy,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:),1);
        a = dat_ens_mean(tind_01:tind_02,:);
        %dat_ukesm(iy,:) = dat_ens_mean(tind_01:tind_02,:);
        dat_ukesm(iy,:) = a(:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        %end
        
        %DJF means
        if iy==1 %skip the first year since don't have December
            dat_annual_DJF(iy,:,:) = NaN*ones([size(dat_ens_mean,2) size(dat_ens_mean,3)]);
            %dat_ukesm_DJF(iy,:,:,:) = NaN;
            dat_annual_ens_DJF(:,iy,:,:) = NaN*ones([size(dat_ens,1) size(dat_ens,3) size(dat_ens,4)]);
        else
            tind_01 = ind_start_offset + (iy-1)*12 + 0; %+1 is Jan, so +0 is Dec
            tind_02 = tind_01+2; %Feb
            dat_annual_DJF(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
            %dat_ukesm_DJF(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
            %calculate the annual means keeping each ensemble separate
            dat_annual_ens_DJF(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        end
        
        % MAM
        tind_01 = ind_start_offset + (iy-1)*12 + 3; %+1 is Jan
        tind_02 = tind_01+2; %May
        dat_annual_MAM(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm_MAM(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens_MAM(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        % JJA
        tind_01 = ind_start_offset + (iy-1)*12 + 6; %+1 is Jan
        tind_02 = tind_01+2; %Aug
        dat_annual_JJA(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm_JJA(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens_JJA(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
        % SON
        tind_01 = ind_start_offset + (iy-1)*12 + 9; %+1 is Jan
        tind_02 = tind_01+2; %Nov
        dat_annual_SON(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm_SON(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens_SON(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
        
    end
    
end


%save(savefile,'/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm.mat','-V7.3'...
%    ,'Nd_annual','Nd_ukesm','Nd_annual_ens','Nd_ens','Nd_ens_mean','Nd_ens_Ndatap','Nd_ens_std','gcm_Plat2D_UM','gcm_Plon2D_UM'...
%    ,'gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM','years_ukesm_1d',);

save(savefile,'-V7.3'); %save all variables


