% For ESGF models from other centres (NCAR-CESM etc.) I think need to run
% this one to read in the NetCDFs and save in the right format :-
% ACSIS_load_ESGF_ensemble_multi_vars_MULTI_model.m

fscale_multi_model_default = 1e-6;

%var_multi_model_DRIVER = 'Nd_clw_weighted_ESGF_no_dz';
var_multi_model_DRIVER = var_ukesm;

had_str = 'HadGEM3_GC31_LL';


switch var_multi_model_DRIVER
    case {'Nd_cf_weighted_UKESM'}
        if iadd_HADGEM==1            
            %for the SW paper rather than the work with Laura
            models = {'HadGEM3-GC31-LL'};            
            %fscale_multi_model_default = 1;
        end
        
    case {'Nd_clw_weighted_ESGF'} 
        if iadd_HADGEM==1            
            %for the SW paper rather than the work with Laura
            models = {'HADGEM3_GC31_LL'};   
            
            had_str = 'HADGEM3_GC31_LL';        
            %fscale_multi_model_default = 1;
        end
    case 'Nd_clw_weighted_ESGF_no_dz_no_ice_total_column_to_zdomain_top'
        %models = {'CESM2',...            
        %    'GFDL-ESM4'};
        
         models = {'BCC-ESM1',...
            'CNRM-CM6-1',...
            'CNRM-ESM2-1',...               
            'MIROC6',...
            'CESM2',...               
            'IPSL-CM6A-LR',...
            'GFDL-CM4',...
            'GFDL-ESM4',...
            'MPI-ESM-1-2-HAM'};
        
        %'EC-Earth3-AerChem',...  
        %'NorESM2-MM',...
        %'NIMS-KMAUKESM1-0-LL',...
        %'MRI-ESM2-0',...
        
    case 'Nd_clw_weighted_ESGF_no_dz'
        
        %Nd models
        models = {'BCC-ESM1',...
            'CNRM-CM6-1',...
            'CNRM-ESM2-1',...
            'EC-Earth3-AerChem',...
            'MPI-ESM-1-2-HAM',...
            'MIROC6',...
            'CESM2',...
            'NorESM2-MM',...
            'NIMS-KMAUKESM1-0-LL',...
            'IPSL-CM6A-LR',...
            'MRI-ESM2-0',...
            'GFDL-ESM4'};
        
        %models = {'IPSL-CM6A-LR'};
        models = {'MPI-ESM-1-2-HAM'};
        
        
        
    case 'reffclwtop'
        
        % var_multi_model_DRIVER = 'reffclwtop';
        % %reff models
        models = {'HadGEM3-GC31-LL',...
            'IPSL-CM6A-LR',...
            'MIROC6',...
            'MIROC-ES2L',...
            'CESM2',...
            'GFDL-ESM4',...
            'NIMS-KMAUKESM1-0-LL',...
            'INM-CM4-8',...
            'INM-CM5-0',...
            'MRI-ESM2-0'...
            'MPI-ESM-1-2-HAM'};
        
        
        
        %models = {'CESM2'};
        
        %'MIROC-ES2L',...
        %'GFDL-CM4',...
        %'GFDL-ESM4'...
        %'CESM2-FV2',...
        %'CESM2-WACCM',...
        %'CESM2-WACCM-FV2',...
        %'MRI-ESM2-0'
        
    case 'cldnvi'
        
        %cldnvi models
        models = {'CESM2-FV2',...
            'CESM2-WACCM',...
            'CESM2-WACCM-FV2',...
            'CESM2',...
            'NIMS-KMAUKESM1-0-LL',...
            'MRI-ESM2-0'...
            'HadGEM3-GC31-LL'};
               
    case {'rsut','rsutcs'}
        if iadd_HADGEM==1            
            %for the SW paper rather than the work with Laura
            models = {'HadGEM3-GC31-LL'};            
            fscale_multi_model_default = 1;
        end        
        
    case 'clt'
       if iadd_HADGEM==1            
            %for the SW paper rather than the work with Laura
            models = {'HadGEM3-GC31-LL'};            
            %fscale_multi_model_default = 1;
       end 
       fscale_multi_model_default = 1e-2;
       
    case {'lwp','lwpic','ts'}
       if iadd_HADGEM==1            
            %for the SW paper rather than the work with Laura
            models = {'HadGEM3-GC31-LL'};            
            %fscale_multi_model_default = 1;
       end 
       fscale_multi_model_default = 1; %1e-2;
       
    case {'od550aer','od550tot'}
       if iadd_HADGEM==1            
            %for the SW paper rather than the work with Laura
            models = {'HADGEM3_GC31_LL'};            
            %fscale_multi_model_default = 1;
            had_str = 'HADGEM3_GC31_LL'; 
       end 
       fscale_multi_model_default = 1; %1e-2;   
       
    otherwise
        fprintf(1,'\nBest to add a case for the variable to deal with e.g. HadGEM3 vs HADGEM3 names, etc');


end


if iload_amip_now==1
    models = {'amip'};
end

for imodel=1:length(models)
    
    
    expt_str = models{imodel}; expt_str2=remove_character(expt_str,'-','_');
    %expt_str = 'multi_model_hist-GHG'; expt_str2='hist_GHG';
    %expt_str = 'multi_model_hist-nat'; expt_str2='hist_nat';
    
    %set defaults
    fscale_multi_model = fscale_multi_model_default;
    var_multi_model = var_multi_model_DRIVER;
    
    switch var_multi_model_DRIVER
        case 'Nd_clw_weighted_ESGF_no_dz_no_ice_total_column_to_zdomain_top'            
            switch expt_str
                %case {'GFDL-ESM4','GFDL-CM4','CESM2','CNRM-CM6-1','CNRM-ESM2-1','IPSL-CM6A-LR','MIROC6','MPI-ESM-1-2-HAM'}
                case {'GFDL-ESM4','GFDL-CM4','CESM2','IPSL-CM6A-LR','MIROC6','MPI-ESM-1-2-HAM'}
                    %The latest no-ice div CF metric with min CF of 0.05.
                    var_multi_model = 'Nd_clw_weighted_ESGF_no_dz_div_CF_no_ice_total_column_to_zdomain_top';                              
            end
            switch expt_str
                case 'BCC-ESM1'
                    fscale_multi_model = 1e6;
                otherwise
                    fscale_multi_model = fscale_multi_model_default;
            end
            
        case 'Nd_clw_weighted_ESGF_no_dz'
            switch expt_str
                case 'BCC-ESM1'
                    fscale_multi_model = 1e6;
                otherwise
                    fscale_multi_model = fscale_multi_model_default;
            end
            
            switch expt_str
                case {'IPSL-CM6A-LR','MIROC6'};%,'MRI-ESM2-0'} %,'MIROC6'
                    %Was on the fence for 'CNRM-CM6-1','CNRM-ESM2-1', but they
                    %look too high when dividing by CF.
                    %'MRI-ESM2-0' didn't get calculated with the normal method for some
                    %reason, but looks too high with div_CF. Can re-do it
                    %without div_cf.
                    %MIROC6 is 3x too high if use div_CF, or 2x too low if
                    %don't...???
                    var_multi_model = 'Nd_clw_weighted_ESGF_no_dz_div_CF';
                %case {'GFDL-ESM4'}
%                    var_multi_model = 'Nd_clw_weighted_ESGF_no_dz_div_CF_total_column_to_zdomain_top';  
                case {'GFDL-ESM4','CESM2'}
                    %The latest no-ice div CF metric with min CF of 0.05.
                    var_multi_model = 'Nd_clw_weighted_ESGF_no_dz_div_CF_no_ice_total_column_to_zdomain_top';
                case {'UKESM1'}
                    var_multi_model = 'Nd_clw_weighted_ESGF_no_dz_no_ice_total_column_to_zdomain_top';
                otherwise
                    var_multi_model = 'Nd_clw_weighted_ESGF_no_dz';
            end
            
        case 'reffclwtop'
            fscale_multi_model = 1e6;
            
        case 'cldnvi'
            fscale_multi_model = 1;
            
    end
    
    % switch var_ukesm
    %     case 'scldncl'
    %         var_multi_model = 'Nd_clw_weighted_ESGF';
    %         fscale_multi_model = 1e-6; %convert frm m^-3 to cm^-3
    %     otherwise
    %         var_multi_model = var_ukesm;
    %         fscale_multi_model = dat_ukesm.fscale;
    % end
    
    %default :-
    load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_' expt_str '_all_' var_multi_model '.mat'];       
    
    switch var_multi_model_DRIVER
        case 'clt'
             switch expt_str
                case 'HadGEM3-GC31-LL' %Special case for HADGEM and rsut that don't have for ESGF data yet.                    
                    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_HADGEM3_GC31_LL_all_clt.mat';                                
             end 
            
        case 'rsut'
            switch expt_str
                case 'HadGEM3-GC31-LL' %Special case for HADGEM and rsut that don't have for ESGF data yet.                                        
                    load_file = '/home/disk/eos15/d.grosvenor/UM/ACSIS_Nd_trends/EGSF_HADGEM-LL_ensemble_timeseries_all_rsut.mat';                  
            end   
            
        case 'rsutcs'
            switch expt_str
                case 'HadGEM3-GC31-LL' %Special case for HADGEM and rsut that don't have for ESGF data yet.                                        
                    load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_HADGEM3_GC31_LL_all_' var_ukesm '.mat'];
            end 
            
            
        case {'Nd_cf_weighted_UKESM','Nd_clw_weighted_UKESM'}
            switch expt_str
                case 'HadGEM3-GC31-LL' %
                    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_HADGEM3_GC31_LL_all_Nd_clw_weighted_ESGF.mat';                
            end  
            
         case 'lwp'
            switch expt_str
                case 'HadGEM3-GC31-LL' %Special case for HADGEM and rsut that don't have for ESGF data yet.                                          
                    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_HADGEM3_GC31_LL_all_lwp.mat';
            end     
            
        case 'lwpic'
             switch expt_str
                case 'HadGEM3-GC31-LL' %Special case for HADGEM and rsut that don't have for ESGF data yet.                                          
                    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_HADGEM3_GC31_LL_all_lwpic.mat';        
             end 
             
        case 'ts'            
            switch expt_str
                case 'HadGEM3-GC31-LL' %Special case for HADGEM and rsut that don't have for ESGF data yet.                                          
                    load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_HADGEM3_GC31_LL_all_' var_ukesm '.mat'];
             end
            
    end
    switch expt_str
        case 'amip'
            load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_AMIP_all_' var_ukesm '.mat'];
                     
    end
    
    load_file_PI = load_file;
    
    %eval_str =['[dat_PI_' expt_str2 ',dat_ukesm_' expt_str2 ',dat_ukesm_' expt_str2 '_DJF,dat_ukesm_' expt_str2 ...
    %    '_JJA,dat_ukesm_' expt_str2 '_MAM,dat_ukesm_' expt_str2 '_SON]='...
    %    'ACSIS_Robson_paper_load_data_generic_FUNC(var_multi_model,expt_str,expt_str2);'];
    eval_str =['[dat_PI_' expt_str2 ',dat_ukesm_' expt_str2 ',dat_ukesm_' expt_str2 '_DJF,dat_ukesm_' expt_str2 ...
        '_JJA,dat_ukesm_' expt_str2 '_MAM,dat_ukesm_' expt_str2 '_SON]='...
        'ACSIS_Robson_paper_load_data_generic_FUNC(load_file,load_file_PI,var_multi_model,fscale_multi_model);'];
    eval(eval_str);
    
    dat_ukesm_str = ['dat_ukesm_' expt_str2];
    eval_str2 = ['[' dat_ukesm_str '.gcm_area_UM] = calc_area_lat_lon2d(' dat_ukesm_str '.gcm_Plat2D_edges_UM,' dat_ukesm_str '.gcm_Plon2D_edges_UM);'];
    eval(eval_str2);
    
    dat = load(load_file,'dat_ens_mean');
    eval_str2 = ['[me_'  expt_str2 ',N_' expt_str2 ',std_' expt_str2 '] = meanNoNan(dat.dat_ens_mean,1);'];
    eval(eval_str2);
    
    %Run script to choose the data for the selected season and then to average
    %over the box region.
    yr_start_trend_box2 = yr_start_trend_box; yr_end_trend_box2 = yr_end_trend_box;
    end_str = ['_' expt_str2];
    ACSIS_Robson_paper_box_means
    
end