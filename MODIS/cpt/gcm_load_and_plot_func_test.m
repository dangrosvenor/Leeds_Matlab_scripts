function [gcm_mplot_stored,implot] = gcm_load_and_plot_func_test(gcm_mplot_stored,implot,am3_dataset2,vars_mplot,days_str_mplot,times_mplot,iload_files,gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly)

    
    switch am3_dataset2
        case 'old'; 
            gcm_savename = 'AM3_old';
        case '2deg'
            gcm_savename = 'AM3_2deg';
        case '0.5deg'
            gcm_savename = 'AM3_0.5deg';
        otherwise
            gcm_savename = am3_dataset2;
    end
    
    if iload_files==1
        ioverride_load_gcm=1;

        load_gcm_process_vars_flag=4;
        %=0 doesn't load in any of them. =1 loads in the 3D files. =2 loads all
        %except 3D files to save memory. =3 is like 2, but doesn't have the height
        %fields and the cloud layer thicknesses to available to load.
        save_gcm_process_vars_flag=0; %also set the save flag

        cosp_flag=0; %this will get overwritten in load_catfiles_timechunks
        cosp_flag4D=0;
        ilwcAPBP=0;  %default of 0 - overwritten in load_catfiles_timechunks

        % ***   select the files (years and GCM type) to load ***

        imod=1;
        %selects files for load_saved_modis_vars.m & modis_make_monthly_averages_multi_year.m
        %

        clear gcm_year_cases gcm_model_cases


%      gcm_year_cases{imod} = 'y1990'; gcm_month_cases{imod} = '01'; gcm_model_cases{imod} = 'AM3'; imod=imod+1; %
%      gcm_year_cases{imod} = 'y1990'; gcm_model_cases{imod} = 'AM3_with_ice'; imod=imod+1; %      
%      gcm_data_case{imod} = 'y2000'; gcm_model_case{imod} = 'AM3';imod=imod+1; %

%       gcm_year_cases{imod} = 'y2001'; gcm_month_cases{imod} = '01'; gcm_model_cases{imod} = 'CAM5'; imod=imod+1; %
%       gcm_year_cases{imod} = 'y2001'; gcm_model_cases{imod} = 'CAM5_with_ice'; imod=imod+1; %       

%      gcm_year_cases{imod} = 'y2006'; gcm_month_cases{imod} = '01'; gcm_model_cases{imod} = 'AM3'; imod=imod+1; %
%      gcm_year_cases{imod} = 'y2006'; gcm_month_cases{imod} = '07'; gcm_model_cases{imod} = 'AM3'; imod=imod+1; %
      
%      gcm_year_cases{imoded} = 'y2007'; gcm_month_cases{imod} = '01'; gcm_model_cases{imod} = 'AM3'; imod=imod+1; %
%      gcm_year_cases{imod} = 'y2007'; gcm_month_cases{imod} = '07'; gcm_model_cases{imod} = 'AM3'; imod=imod+1; %      
      


filedir_gcm_load = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/';

%new style of chunkded data (data in time chunks)
% ************************************************
load_catfiles_timechunks
% ************************************************
load_gcm_processed_data



    end
    
    
    
    
%%now plot and save

for ivar_mplot = 1:length(vars_mplot)
    
    for idays_mplot = 1:length(days_str_mplot)

    for itime_mplot = 1:length(times_mplot)
        implot=implot+1; %overall counter

        mod_data_type='GCM';   gcm_time_of_day_select=2;
        modis_data_plot= vars_mplot{ivar_mplot};
        ioverride_plotglobal_thresh=1;

        time_mean_str = days_str_mplot{idays_mplot};
        % --------------------
           days_for_season;  %Selects the required days of the year given time_mean_str
        % --------------------        
        times_required = times_mplot{itime_mplot};
        ioverride_time_selection=1;
        
 
%        proj_type='polar';
        proj_type='global oval'; %select this for anything but polar stereographic
        data_select='specific_modis_data';
        ifilter_ndays=1 %flag for whether to cut out points for which there aren't many days
        icontour=0;
        %time_mean_str = 'ALL';        
        noplot=0;
        irestrict_domain=1 %whether to restrict the domain or not (use for regional domains)
        
        plot_global_maps
        if exist('mean_bias_gcm')
            gcm_mplot_stored(implot).mean_bias = mean_bias_gcm;
            gcm_mplot_stored(implot).RMSE = RMSE_gcm;
            gcm_mplot_stored(implot).gcm_str = gcm_str;
            gcm_mplot_stored(implot).days = days_str_mplot{idays_mplot};
            gcm_mplot_stored(implot).times = times_mplot{itime_mplot};
        end
        
        
        saveas_ps_fig_emf(gcf,[savedir_mplot modis_data_plot '_' gcm_savename '_' time_mean_str],'',0,1);


    end
    
    end

end

    
    