%To save data from these loops might want to use
%   save_vars_mat_run.m

try
    
if ~exist('ioverride_DRIVER') | ioverride_DRIVER==0
    fid_mtab = fopen('~/MOD_CAL_mbias_table.txt','wt');
    fid_Rtab = fopen('~/MOD_CAL_Rbias_table.txt','wt');
end

%Some defaults that need setting in plot_global_maps
plot_global_maps_defaults

%save_file='~/MODIS_vs_AMSRE_LWP_bias_saved.mat';
%save_file='~/MODIS_effect_of_CF_screening_Nd_saved.mat';

%fid_mean_tab = 
isave_vals=0; %flag to say whether we want to save some values to a file
isave_plot=1; %whether to save the plot
savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/';

proj_type='global oval'; %select this for anything but polar stereographic
data_select='specific_modis_data';
ifilter_ndays=0; %flag for whether to cut out points for which there aren't many days
icontour=0;
noplot=0;
plot_num_datapoints=0;
irestrict_domain = 1;
thresh_LAT = [-40 10];  thresh_LON = [-140 -50]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
%thresh_LAT = [40 70];  thresh_LON = [-40 10]; lon_ticks=[-40:5:10]; lat_ticks=[40:10:90]; %Iceland volcano (Jan 2015)
%thresh_LAT = [40 70];  thresh_LON = [-40 -12]; lon_ticks=[-40:5:10]; lat_ticks=[40:10:90]; %Iceland volcano - Atlantic only region (Jan 2015)
%thresh_LAT = [40 89];  thresh_LON = [-40 10]; lon_ticks=[-40:5:10]; lat_ticks=[40:10:90]; %Iceland volcano (Jan 2015)
%thresh_LAT = [40 90];  thresh_LON = [-180 180]; lon_ticks=[-40:5:10]; lat_ticks=[40:10:90]; %Iceland volcano (Jan 2015)

thresh_zeroCF = 0.05;

plot_set = 'MODIS vs CALIPSO bias';
%plot_set = 'MODIS vs AMSRE LWP';
%plot_set = 'MODIS vs AMSRE LWP bias';
%plot_set = 'MODIS vs AMSRE LWP bias with Ndays and mean CF contours';
% plot_set = 'MODIS Nd - effect of screening';
% plot_set = 'Model minus MODIS Nd';
% %plot_set = 'Model Nd daytime minus nighttime';
% plot_set = 'MODIS Nd multiple year seasonal cycles';
%plot_set = 'MODIS Nd one year all wavelengths seasonal cycles';
%plot_set = 'MODIS Nd one year all wavelengths seasonal cycles from individual swaths';
%plot_set = 'MODIS Nd diffs between datasets one year all wavelengths seasonal cycles from individual swaths';
    %E.g. difference between the mockL3 and actual L3 as looked into for
    %Dan McCoy in May2014.
%plot_set = 'model_LWP_vs_AMSRE';
%plot_set = 'model_LWP_vs_AMSRE_conditional_sampling';  %can also switch off conditional sampling here.
%plot_set = 'model_TLWP_vs_AMSRE_conditional_sampling';
%plot_set = 'model_TLWP_PDF'; %Also works for LWP - select which within the file
%plot_set = 'model_precip_rate_PDF'; %Using this to calcualte the dirunal precip range for the scatter plot as of 25th Aug, 2015
%plot_set = 'DRIVER_UM_plots';
%plot_set = 'multi_plot_many_years_for_one_season';  %(From MODIS L3) - map plots of e.g. Nd
%plot_set = 'multi_plot_many_years_for_one_season_reff';  %(From MODIS L3) - map plots of e.g. Nd
%plot_set = 'multi_plot_many_years_for_one_season_Nd';  %(From MODIS L3) - map plots of e.g. Nd
%plot_set = 'multi_plot_many_years_for_one_season_Nd_UKESM';  %(From MODIS L3) - map plots of e.g. Nd
%plot_set = 'multi_plot_many_years_for_one_season_CTT_PDFs';  %(From MODIS L3)
%plot_set = 'multi_plot_many_years_for_one_season_reff_PDFs';  %(From MODIS L3)
%plot_set = 'multi_plot_many_years_for_one_season_reff_anomalies';  %(From MODIS L3) - map plots of e.g. Nd
%plot_set = 'Florent Iceland Nd data for Aqua QA Jan 2016';  % Aqua QA Nd dataset for Iceland volcano paper (previously only
   % did this for Terra). For Jan 2016 revisions of Haywood paper.


switch plot_set
    case 'multi_plot_many_years_for_one_season_Nd_UKESM'
        multi_plot_many_years_for_one_season_Nd2_QA_ndays  
    case 'Florent Iceland Nd data for Aqua QA Jan 2016'
        multi_plot_many_years_one_season_Nd2_QA_ndays_FloJan2016Iceland
    case 'multi_plot_many_years_for_one_season_Nd'
        multi_plot_many_years_for_one_season_Nd2        
    case 'multi_plot_many_years_for_one_season_CTT_PDFs'
        multi_plot_many_years_for_one_season_CTT_PDFs
 case 'multi_plot_many_years_for_one_season_reff_PDFs'
        multi_plot_many_years_for_one_season_reff_PDFs        
     case 'multi_plot_many_years_for_one_season_reff'
        multi_plot_many_years_for_one_season_reff2       
    case 'multi_plot_many_years_for_one_season'
        multi_plot_many_years_for_one_season
    case 'DRIVER_UM_plots'
        DRIVER_UM_plots
    case 'model_TLWP_PDF'
        model_TLWP_PDF
    case 'model_precip_rate_PDF'
        model_precip_rate_PDF        
    case 'model_TLWP_vs_AMSRE_conditional_sampling'
        model_TLWP_vs_AMSRE_conditional_sampling
    case 'model_LWP_vs_AMSRE_conditional_sampling'
        model_LWP_vs_AMSRE_conditional_sampling
    case 'model_LWP_vs_AMSRE'
        model_LWP_vs_AMSRE
    case 'MODIS Nd vs height various months from individual swaths CF80 no_seaice'
        MODIS_Nd_vs_height_individ_swaths_CF80_no_seaice
    case 'MODIS Nd one year all wavelengths seasonal cycles from individual swaths'
        MODIS_Nd1Y_AllWlengths_seasonal_individ_swaths
        
    case 'MODIS Nd diffs between datasets one year all wavelengths seasonal cycles from individual swaths'
        MODIS_Nd_diffs_between_datasets_seasonal_cycles_L3
        
    case 'MODIS Nd one year all wavelengths seasonal cycles'
       % Plot seasonal cycles of MODIS Nd separately for different wavelengths using 2007 mock L3 data.
       % Also split up in terms of CF.
       % Saved this data using save_vars_mat_run
       
       fraction_days_plot = 0; %flag for wheterh to do the fraction of days remaining after screening plot
       
%        modis_data_plots={'MODIS LWP minus AMSRE'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';
%        modis_data_plots={'Number of droplets cell values time mean - specific days'}; fraction_plotcase='Fraction of days remaining after screening';
%         modis_data_plots={'LWP AMSRE time3'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';

       % -- Specify the plots to call in plot_global etc.
%        modis_data_plots={'Dummy string'}; %just put somethign in here to make the loop work
         modis_data_plots ={'Number of droplets 1.6um cell values time mean - specific days', ...
             'Number of droplets cell values time mean - specific days', ...
             'Number of droplets 3.7um cell values time mean - specific days', ...
             };
    
%        zero_CF_out_of_range=0;
        time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';
        
        set_screening = {'none','NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
            %Mean CTH:-
        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
            %Max CTH:-
%        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};        
            %Mean CTH with sea-ice screening:-
        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
%        set_screening = {'Dummy string'}; %Dummy string to enter loop, but will set this later - acutally might not 
        %be used for the GCM plot
        
        set_screening = {'none'}; %if using mockL3 would generally use this as the screening has already
           %been done. However, if using concatenated L3 from separate overpasses then would want 
           %to set some screening


%set a default screening - then can change individual ones in the loop
% screening_eval_str=[ ...
%         'thresh_ndays=1;' ...
%         'thresh_SZA=[0 65];' ...
%         'thresh_CF=[0.8 1.00001];' ...
%         'thresh_NP=50;' ...
%         'thresh_sensZA=[0 90];'     ...
%         'thresh_CTT = [273 273+100];'  ...
%         'thresh_relAZ = [0 180];'   ...      
%         'thresh_CTH = [-0.01 3.2];' ...
%         'thresh_stdW = [0 1e9];' ...
%         'thresh_sigCTT = [0 1e9];' ...
%         'thresh_reff = [0 14];' ...
%         ];
    
    screening_eval_str='';

    %just ones that alter it from the above baseline
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];'};
screen_edits = {'thresh_CF=[0.8 1.00001]','thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];','thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];','thresh_CF=[0.8 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_reff=[0 30];'};
screen_edits = {''};

screen_edits = {'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %all CF, all re
%screen_edits = {'thresh_reff=[0 30];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %CF>80 and all CF (all re)

 screen_edits = {'thresh_CF=[0.1 0.4]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
     ,'thresh_CF=[0.4 0.6]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
     ,'thresh_CF=[0.6 0.8]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
     ,'thresh_CF=[0.8 1.00001]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...    
     };  %ranges of CF (applies to model too when choosing re CF range)
 
 screen_edits = {''};
 
 set_years = [2007];

        ifilter_ndays=0;
        
       
        
        
        set_MOD = {'MOD06 liquid','MOD35'};
         set_MOD = {'MOD35'};  %redundant for Nd
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
        set_real_abs = {'absolute','percentage'};
        set_real_abs = {'percentage'};   %redundant for Nd
%        set_real_abs = {'absolute'};
        
        %what to plot for the black contours. - %redundant for Nd
        cont_dat_choose = 'calipso mid+highCF';  
        
        %what to use for filtering - %redundant for Nd
        filt_dat_choose = 'calipso mid+highCF';
        ifilter_clhcalipso=1; %whether to do the filtering or not
        thresh_clh = [-0.01 0.3];
        
        iuse_MODIS_CF_error=0; %if =1 then uses both the AMSRE comparison error and the change due to changing CF
        %from MODIS to plot the contours for the more dubious area.
        %Otherwise just uses the LWP error. Decided that the latter is
        %better as there will be a real component ot the change in Nd with
        %CF.
        
        iuse_saved_mask=0; %redundant for Nd
        saved_mask_name = 'imask_ndays15_cal_midhigh';
        
        clear RMSE_DRIVER mean_bias_DRIVER season_mean_driver season_Ndatap_driver season_std_driver mean_sig2_driver Nspatial_driver season_timestd_driver season_datatype_driver season_vals_timelabs_driver screen_type_driver
                        
        icount_DRIVER=0;             
        
        for iloop_plot=1:length(modis_data_plots)
            modis_data_plot = modis_data_plots{iloop_plot};

            for iscreen_set=1:length(set_screening)
                screen_type = set_screening{iscreen_set};

                for iscreen_thresh=1:length(screen_edits)
                    eval(screening_eval_str); %eval the basline screening
                    eval(screen_edits{iscreen_thresh}); %and now the edits to the baseline

                    for iloop_MOD=1:length(set_MOD)
                        MODCF_type = set_MOD{iloop_MOD}

                        for iloop_abs=1:length(set_real_abs)
                            error_type= set_real_abs{iloop_abs};
                                for iloop_years=1:length(set_years)
                                    years_required_for_mean = set_years(iloop_years)
                                    icount_DRIVER = icount_DRIVER+1;
                            

            %             modis_data_plot = modis_data_plots{};  %set above



                            ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
%                            modis_data_plot = modis_data_plots{iloop_plot};

                            ioverride_years_time_screen=1; %required for the different years
                            ioverride_monthly_options=1;
                                  multi_case = 'Monthly';
                                  plot_script = 'plot_global_maps';
                            ioverride_location_selection2=1;
                                  LAT_vals = [-60 -40;]; LON_vals = [-180 179;]; %Southern Ocean
                            

                                    iset_min_clim=0;
                                    clim_min=0;
                                    iset_max_clim=0;
                                    %        clim_max=225;
                                    clim_max=150;
                                    
%                            times_required = [0:24]; %need to specify all hours (in 1 hour increments)
%                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
     % The above should be set in monthly_means_from_plot_global                            
                            


                            
                                                        
                           % plot_global_maps
                           monthly_means_from_plot_global
%                            Psave_MODIS = P;
%                            Npoints_save_MODIS = Npoints;
%                            thresh_str_MODIS = thresh_str;
%                            mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
%                            RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;

                          thresh_str_driver{icount_DRIVER} = thresh_str; %save the thresh_str used
                          screen_type_driver{icount_DRIVER} = screen_type;
                          years_required_for_mean_driver{icount_DRIVER} = years_required_for_mean;

                        season_mean_driver{icount_DRIVER} = season_mean;
                        season_Ndatap_driver{icount_DRIVER} = season_Ndatap;
                        season_std_driver{icount_DRIVER} = season_std;
                        mean_sig2_driver{icount_DRIVER} = mean_sig2;
                        Nspatial_driver{icount_DRIVER} = Nspatial;
                        season_timestd_driver{icount_DRIVER} = season_timestd;

                        season_datatype_driver{icount_DRIVER} = season_datatype;
                        season_vals_timelabs_driver{icount_DRIVER} = season_vals_timelabs;
                        
                        if exist('thresh_str_multiL2')
                            thresh_str_multiL2_driver{icount_DRIVER} = thresh_str_multiL2;
                        else
                            thresh_str_multiL2_driver{icount_DRIVER} = '';
                        end

                        if exist('LON_str')
                            season_LON_str{iseason} = LON_str;
                            season_LAT_str{iseason} = LAT_str;
                        end
                        
                        
                            
%                            if isave_vals==1
%                                fprintf(fid_mtab,'%f ',mean_bias_MODIS);
%                               fprintf(fid_Rtab,'%f ',RMSE_MODIS);
%                            end
                          
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
                            
                        end

                    end

                end

            end

            end
        end
        
        
        
        case 'MODIS Nd multiple year seasonal cycles'        
       % Plot seasonal cycles of MODIS Nd separately for different years.
       % Also split up in terms of CF.
       
       fraction_days_plot = 0; %flag for wheterh to do the fraction of days remaining after screening plot
       
%        modis_data_plots={'MODIS LWP minus AMSRE'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';
%        modis_data_plots={'Number of droplets cell values time mean - specific days'}; fraction_plotcase='Fraction of days remaining after screening';
%         modis_data_plots={'LWP AMSRE time3'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';

        modis_data_plots={'Dummy string'}; %just put somethign in here to make the loop work
    
%        zero_CF_out_of_range=0;
        time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';
        
        set_screening = {'none','NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
%        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};        
%        set_screening = {'Dummy string'}; %Dummy string to enter loop, but will set this later - acutally might not 
        %be used for the GCM plot

%set a default screening - then can change individual ones in the loop
screening_eval_str=[ ...
        'thresh_ndays=1;' ...
        'thresh_SZA=[0 65];' ...
        'thresh_CF=[0.8 1.00001];' ...
        'thresh_NP=50;' ...
        'thresh_sensZA=[0 90];'     ...
        'thresh_CTT = [273 273+100];'  ...
        'thresh_relAZ = [0 180];'   ...      
        'thresh_CTH = [-0.01 3.2];' ...
        'thresh_stdW = [0 1e9];' ...
        'thresh_sigCTT = [0 1e9];' ...
        'thresh_reff = [0 14];' ...
        ];

    %just ones that alter it from the above baseline
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];'};
screen_edits = {'thresh_CF=[0.8 1.00001]','thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];','thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];','thresh_CF=[0.8 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_reff=[0 30];'};
screen_edits = {''};

screen_edits = {'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %all CF, all re
%screen_edits = {'thresh_reff=[0 30];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %CF>80 and all CF (all re)

 screen_edits = {'thresh_CF=[0.1 0.4]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
     ,'thresh_CF=[0.4 0.6]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
     ,'thresh_CF=[0.6 0.8]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...
     ,'thresh_CF=[0.8 1.00001]; thresh_reff=[0 30]; thresh_CTT = [263 273+100];'...    
     };  %ranges of CF (applies to model too when choosing re CF range)
 
 set_years = [2006:2010];

        ifilter_ndays=1;
        
       
        
        
        set_MOD = {'MOD06 liquid','MOD35'};
         set_MOD = {'MOD35'};  %redundant for Nd
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
        set_real_abs = {'absolute','percentage'};
        set_real_abs = {'percentage'};   %redundant for Nd
%        set_real_abs = {'absolute'};
        
        %what to plot for the black contours. - %redundant for Nd
        cont_dat_choose = 'calipso mid+highCF';  
        
        %what to use for filtering - %redundant for Nd
        filt_dat_choose = 'calipso mid+highCF';
        ifilter_clhcalipso=1; %whether to do the filtering or not
        thresh_clh = [-0.01 0.3];
        
        iuse_MODIS_CF_error=0; %if =1 then uses both the AMSRE comparison error and the change due to changing CF
        %from MODIS to plot the contours for the more dubious area.
        %Otherwise just uses the LWP error. Decided that the latter is
        %better as there will be a real component ot the change in Nd with
        %CF.
        
        iuse_saved_mask=0; %redundant for Nd
        saved_mask_name = 'imask_ndays15_cal_midhigh';
        
        clear RMSE_DRIVER mean_bias_DRIVER season_mean_driver season_Ndatap_driver season_std_driver mean_sig2_driver Nspatial_driver season_timestd_driver season_datatype_driver season_vals_timelabs_driver
                        
        icount_DRIVER=0;             
        
        for iloop_plot=1:length(modis_data_plots)
            modis_data_plot = modis_data_plots{iloop_plot};

            for iscreen_set=1:length(set_screening)
                screen_type = set_screening{iscreen_set};

                for iscreen_thresh=1:length(screen_edits)
                    eval(screening_eval_str); %eval the basline screening
                    eval(screen_edits{iscreen_thresh}); %and now the edits to the baseline

                    for iloop_MOD=1:length(set_MOD)
                        MODCF_type = set_MOD{iloop_MOD}

                        for iloop_abs=1:length(set_real_abs)
                            error_type= set_real_abs{iloop_abs};
                                for iloop_years=1:length(set_years)
                                    years_required_for_mean = set_years(iloop_years)
                                    icount_DRIVER = icount_DRIVER+1;
                            
                    % -- run MODIS Nd first        
                            modis_data_plot ='Number of droplets cell values time mean - specific days';
%                            screen_type = 'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff';


                            ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
%                            modis_data_plot = modis_data_plots{iloop_plot};

                            ioverride_years_time_screen=1; %required for the different years
                            ioverride_monthly_options=1;
                                  multi_case = 'Monthly';
                                  plot_script = 'plot_global_maps';
                            ioverride_location_selection2=1;
                                  LAT_vals = [-60 -40;]; LON_vals = [-180 169;]; %Southern Ocean
                            

                                    iset_min_clim=0;
                                    clim_min=0;
                                    iset_max_clim=0;
                                    %        clim_max=225;
                                    clim_max=150;
                                    
%                            times_required = [0:24]; %need to specify all hours (in 1 hour increments)
%                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
     % The above should be set in monthly_means_from_plot_global                            
                            


                            
                                                        
                           % plot_global_maps
                           monthly_means_from_plot_global
%                            Psave_MODIS = P;
%                            Npoints_save_MODIS = Npoints;
%                            thresh_str_MODIS = thresh_str;
%                            mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
%                            RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;

                          thresh_str_driver{icount_DRIVER} = thresh_str; %save the thresh_str used
                          years_required_for_mean{icount_DRIVER} = years_required_for_mean;

                        season_mean_driver{icount_DRIVER} = season_mean;
                        season_Ndatap_driver{icount_DRIVER} = season_Ndatap;
                        season_std_driver{icount_DRIVER} = season_std;
                        mean_sig2_driver{icount_DRIVER} = mean_sig2;
                        Nspatial_driver{icount_DRIVER} = Nspatial;
                        season_timestd_driver{icount_DRIVER} = season_timestd;

                        season_datatype_driver{icount_DRIVER} = season_datatype;
                        season_vals_timelabs_driver{icount_DRIVER} = season_vals_timelabs;

                        if exist('LON_str')
                            season_LON_str{iseason} = LON_str;
                            season_LAT_str{iseason} = LAT_str;
                        end
                        
                        
                            
%                            if isave_vals==1
%                                fprintf(fid_mtab,'%f ',mean_bias_MODIS);
%                               fprintf(fid_Rtab,'%f ',RMSE_MODIS);
%                            end
                          
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
                            
                        end

                    end

                end

            end

            end
        end
        
        
    case 'Model minus MODIS Nd'        
       % Creates maps of mean difference in Nd between model and MODIS
       % (based on the time-means for each location).
       % Taking the approach of running plot_global once for MODIS and then once
       % for the model since they use different screenings. Then plotting
       % the difference.
       
       fraction_days_plot = 0; %flag for wheterh to do the fraction of days remaining after screening plot
       
%        modis_data_plots={'MODIS LWP minus AMSRE'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';
%        modis_data_plots={'Number of droplets cell values time mean - specific days'}; fraction_plotcase='Fraction of days remaining after screening';
%         modis_data_plots={'LWP AMSRE time3'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';

        modis_data_plots={'Dummy string'}; %just put somethign in here to make the loop work
    
%        zero_CF_out_of_range=0;
        time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';
        
        set_screening = {'none','NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
%        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};        
%        set_screening = {'Dummy string'}; %Dummy string to enter loop, but will set this later - acutally might not 
        %be used for the GCM plot

%set a default screening - then can change individual ones in the loop
screening_eval_str=[ ...
        'thresh_ndays=15;' ...
        'thresh_SZA=[0 90];' ...
        'thresh_CF=[0.8 1.00001];' ...
        'thresh_NP=50;' ...
        'thresh_sensZA=[0 90];'     ...
        'thresh_CTT = [273 273+100];'  ...
        'thresh_relAZ = [0 180];'   ...      
        'thresh_CTH = [-0.01 3.2];' ...
        'thresh_stdW = [0 1e9];' ...
        'thresh_sigCTT = [0 1e9];' ...
        'thresh_reff = [0 14];' ...
        ];

    %just ones that alter it from the above baseline
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];'};
screen_edits = {'thresh_CF=[0.8 1.00001]','thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];','thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];','thresh_CF=[0.8 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_reff=[0 30];'};
screen_edits = {''};

screen_edits = {'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %all CF, all re
%screen_edits = {'thresh_reff=[0 30];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};  %CF>80 and all CF (all re)
% screen_edits = {'thresh_CF=[0.1 0.4]; thresh_reff=[0 30];'...
%     ,'thresh_CF=[0.4 0.6]; thresh_reff=[0 30];'...
%     ,'thresh_CF=[0.6 0.8]; thresh_reff=[0 30];'...
%     ,'thresh_CF=[0.8 1.00001]; thresh_reff=[0 30];'...    
%     };  %ranges of CF (applies to model too when choosing re CF range)

        ifilter_ndays=1;
        
       
        
        
        set_MOD = {'MOD06 liquid','MOD35'};
         set_MOD = {'MOD35'};  %redundant for Nd
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
        set_real_abs = {'absolute','percentage'};
        set_real_abs = {'percentage'};   %redundant for Nd
%        set_real_abs = {'absolute'};
        
        %what to plot for the black contours. - %redundant for Nd
        cont_dat_choose = 'calipso mid+highCF';  
        
        %what to use for filtering - %redundant for Nd
        filt_dat_choose = 'calipso mid+highCF';
        ifilter_clhcalipso=1; %whether to do the filtering or not
        thresh_clh = [-0.01 0.3];
        
        iuse_MODIS_CF_error=0; %if =1 then uses both the AMSRE comparison error and the change due to changing CF
        %from MODIS to plot the contours for the more dubious area.
        %Otherwise just uses the LWP error. Decided that the latter is
        %better as there will be a real component ot the change in Nd with
        %CF.
        
        iuse_saved_mask=0; %redundant for Nd
        saved_mask_name = 'imask_ndays15_cal_midhigh';
        
        clear RMSE_DRIVER mean_bias_DRIVER
        icount_DRIVER=0;             
        
        for iloop_plot=1:length(modis_data_plots)
            modis_data_plot = modis_data_plots{iloop_plot};

            for iscreen_set=1:length(set_screening)
                screen_type = set_screening{iscreen_set};

                for iscreen_thresh=1:length(screen_edits)
                    eval(screening_eval_str); %eval the basline screening
                    eval(screen_edits{iscreen_thresh}); %and now the edits to the baseline

                    for iloop_MOD=1:length(set_MOD)
                        MODCF_type = set_MOD{iloop_MOD}

                        for iloop_abs=1:length(set_real_abs)
                            error_type= set_real_abs{iloop_abs};
                            icount_DRIVER = icount_DRIVER+1;
                            
                    % -- run MODIS Nd first        
                            modis_data_plot ='Number of droplets cell values time mean - specific days';
%                            screen_type = 'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff';


                            ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
%                            modis_data_plot = modis_data_plots{iloop_plot};
                            

                                    iset_min_clim=0;
                                    clim_min=0;
                                    iset_max_clim=0;
                                    %        clim_max=225;
                                    clim_max=150;
                                    
                            times_required = [0:24]; %need to specify all hours (in 1 hour increments)
                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                            
                                                        
                            plot_global_maps
                            Psave_MODIS = P;
                            Npoints_save_MODIS = Npoints;
                            thresh_str_MODIS = thresh_str;
%                            mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
%                            RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;
                            
%                            if isave_vals==1
%                                fprintf(fid_mtab,'%f ',mean_bias_MODIS);
%                               fprintf(fid_Rtab,'%f ',RMSE_MODIS);
%                            end
                          
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
                            
% --
                    % -- now do the model Nd      
                            modis_data_plot ='Max Nd in all non-screened layers GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
%                            modis_data_plot ='Max Nd noCF GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
%                            modis_data_plot ='Max Nd in low cloud with cloud screening'; mod_data_type='GCM';  gcm_time_of_day_select=2;                                                        
%                            modis_data_plot ='Max Nd in low cloud without cloud screening'; mod_data_type='GCM';  gcm_time_of_day_select=2;                            
                            
                            %tells 'Max Nd in low cloud without cloud screening' what screening to do.
                            screen_using_cosp = 're CF range'; 
                            screen_using_cosp = 're';
                            screen_using_cosp = 'none';                            


                            ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            
                           
                            % -- N.B. - If using inew_cticks the
                            % differencing later on won't work (since the P
                            % values are scaled). Maybe save the originals
                            % before scaling and use them?
                                    iset_min_clim=0;
                                    clim_min=0;
                                    iset_max_clim=0;
                                    %        clim_max=225;
                                    clim_max=300;                              
                            
                                                        
                            plot_global_maps
                            %interpolate the model field to the MODIS grid
                            %Uses P_save because this is the P field before
                            %any scaling for the colobar etc.
                            P_int =  griddata(Plon2D,Plat2D,P_save,MLON,MLAT');                            
                            Psave_model = P_int;
                            Npoints_save_model = Npoints;
%                            mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
%                            RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;

                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
                            
%  ---   Plot the number of days with data for the model (to
%                      gauge how important a certain cloud fraction is)

                            colormap_choose=jet; %default

                            P = Ndays2./size(dat_modis,1); %number of days used out the total possible
%                            Npoints = min(Npoints_save_model,Npoints_save_MODIS);

                            MODIS_varname2_plot = ['Fraction of days with model N_d for' cosp_screen_str];

                            units_str_plot='';
                            title_info = '';
                            
                            %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=1;
                            %load in the bias data
                            save_file='~/MODIS_vs_AMSRE_LWP_bias_saved.mat';
                            load(save_file,'Psave'); %don't want to load screen_edits
                            P_LWP_bias_allCF_reLT14 = Psave{4}; %Precentage bias
                            

                            
                            if iuse_MODIS_CF_error==1

                                save_file='~/MODIS_effect_of_CF_screening_Nd_saved.mat';
                                load(save_file,'Psave'); %don't want to load screen_edits
                                P_Nd_bias_allCF_vs_highCF = 100*(Psave{1}-Psave{2})./Psave{2};
                                %Psave{1} is CF>0.8 and Psave{2} is all CF.

                                %Choose the max out of the two percentage
                                %values.
                                cont_dat = max(abs(P_LWP_bias_allCF_reLT14),abs(P_Nd_bias_allCF_vs_highCF));

                            else
                                cont_dat = abs(P_LWP_bias_allCF_reLT14);
                            end
                            %interpolate the error data to match the GCM grid
                            cont_dat =  griddata(MLON,MLAT',cont_dat,Plon2D,Plat2D);  
                            
                            %Plot contour of thresh_P percent
                            thresh_P = 20;
                            cont_ints=[thresh_P thresh_P];
                            
                            inew_cticks=0;


                            iset_min_clim=1;
                            clim_min=0;
                            iset_max_clim=1;
                            clim_max=0.5;

                            ifilter_ndays=0;
                            thresh_ndays=15;


                            ioverride_plotglobal_thresh=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
%                            mod_data_type='timeseries3 lambert';
                            modis_data_plot = 'Generic plot specified outside of script';
                            mod_data_type='GCM';  gcm_time_of_day_select=2; 

                            plot_global_maps
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end

                            

                            
                            
                            
%  ---               Now plot the difference

                            colormap_choose=jet; %default

                            P = 100*(Psave_model-Psave_MODIS)./Psave_MODIS;
                            Npoints = min(Npoints_save_model,Npoints_save_MODIS);

                            MODIS_varname2_plot = ['Percentage diff in N_d for model (' gcm_str ' using ' modis_data_plot cosp_screen_str ') vs MODIS for ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str_MODIS];

                            units_str_plot='%';
                            title_info = '';
                            
                            %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=1;
                            %load in the bias data
                            save_file='~/MODIS_vs_AMSRE_LWP_bias_saved.mat';
                            load(save_file,'Psave'); %don't want to load screen_edits
                            P_LWP_bias_allCF_reLT14 = Psave{4}; %Precentage bias
                            

                            
                            if iuse_MODIS_CF_error==1

                                save_file='~/MODIS_effect_of_CF_screening_Nd_saved.mat';
                                load(save_file,'Psave'); %don't want to load screen_edits
                                P_Nd_bias_allCF_vs_highCF = 100*(Psave{1}-Psave{2})./Psave{2};
                                %Psave{1} is CF>0.8 and Psave{2} is all CF.

                                %Choose the max out of the two percentage
                                %values.
                                cont_dat = max(abs(P_LWP_bias_allCF_reLT14),abs(P_Nd_bias_allCF_vs_highCF));

                            else
                                cont_dat = abs(P_LWP_bias_allCF_reLT14);
                            end
                            %Plot contour of thresh_P percent
                            thresh_P = 20;
                            cont_ints=[thresh_P thresh_P];
                            
                            inew_cticks=0;


                            iset_min_clim=1;
                            clim_min=-50;
                            iset_max_clim=1;
                            clim_max=50;

                            ifilter_ndays=0;
                            thresh_ndays=15;


                            ioverride_plotglobal_thresh=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
                            modis_data_plot = 'Generic plot specified outside of script';

                            plot_global_maps
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end



                        end

                    end

                end

            end

        end
           
        if isave_vals==1
            fprintf(fid_mtab,'\n');
            fprintf(fid_Rtab,'\n');
        end

     
if ~exist('ioverride_DRIVER') | ioverride_DRIVER==0        
        if isave_vals==1
            fclose(fid_mtab);
            fclose(fid_Rtab);
        end
end
        
        

        
                          

     case 'MODIS Nd - effect of screening'    
       % Creates maps of mean Nd for no CF screening and for CF>0.8 only
       % and computes the difference
       % Need to have two screenings running in the loop (as does the
       % difference)
       
       icontour_set=0; %wheter to plot contours or not
       
       fraction_days_plot = 0; %flag for wheterh to do the fraction of days remaining after screening plot
       
%        modis_data_plots={'MODIS LWP minus AMSRE'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';
        modis_data_plots={'Number of droplets cell values time mean - specific days'}; fraction_plotcase='Fraction of days remaining after screening';
%         modis_data_plots={'LWP AMSRE time3'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';

%        zero_CF_out_of_range=0;
        time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';
        
        set_screening = {'none','NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};


%set a default screening - then can change individual ones in the loop
screening_eval_str=[ ...
        'thresh_ndays=15;' ...
        'thresh_SZA=[0 90];' ...
        'thresh_CF=[0.8 1.00001];' ...
        'thresh_NP=50;' ...
        'thresh_sensZA=[0 90];'     ...
        'thresh_CTT = [273 273+100];'  ...
        'thresh_relAZ = [0 180];'   ...      
        'thresh_CTH = [-0.01 3.2];' ...
        'thresh_stdW = [0 1e9];' ...
        'thresh_sigCTT = [0 1e9];' ...
        'thresh_reff = [0 14];' ...
        ];

    %just ones that alter it from the above baseline
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];'};
screen_edits = {'thresh_CF=[0.8 1.00001]','thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];','thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];','thresh_CF=[0.8 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_reff=[0 30];'};
screen_edits = {''};
screen_edits = {'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
%screen_edits = {'thresh_reff=[0 30];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
%screen_edits = {'thresh_reff=[0 30]; thresh_CF=[-0.01 0.2];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
screen_edits = {'thresh_reff=[0 30]; thresh_CF=[0.2 0.4];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
%screen_edits = {'thresh_reff=[0 30]; thresh_CF=[0.4 0.6];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
%screen_edits = {'thresh_reff=[0 30]; thresh_CF=[0.6 0.8];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
%screen_edits = {'thresh_reff=[0 30]; thresh_CF=[0.8 1.01];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
%screen_edits = {'thresh_reff=[0 30]; thresh_CF=[0.2 0.4];', 'thresh_CF=[0.8 1.00001]; thresh_reff=[0 30];'};

        ifilter_ndays=1;
        
       
        
        
        set_MOD = {'MOD06 liquid','MOD35'};
         set_MOD = {'MOD35'};  %redundant for Nd
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
        set_real_abs = {'absolute','percentage'};
        set_real_abs = {'percentage'};   %redundant for Nd
%        set_real_abs = {'absolute'};
        
        %what to plot for the black contours. - %redundant for Nd
        cont_dat_choose = 'calipso mid+highCF';  
        
        %what to use for filtering - %redundant for Nd
        filt_dat_choose = 'calipso mid+highCF';
        ifilter_clhcalipso=1; %whether to do the filtering or not
        thresh_clh = [-0.01 0.3];
        
        iuse_saved_mask=0; %redundant for Nd
        saved_mask_name = 'imask_ndays15_cal_midhigh';
        
        clear RMSE_DRIVER mean_bias_DRIVER Psave
        icount_DRIVER=0;             
        
        for iloop_plot=1:length(modis_data_plots)
            modis_data_plot = modis_data_plots{iloop_plot};

            for iscreen_set=1:length(set_screening)
                screen_type = set_screening{iscreen_set};

                for iscreen_thresh=1:length(screen_edits)
                    eval(screening_eval_str); %eval the basline screening
                    eval(screen_edits{iscreen_thresh}); %and now the edits to the baseline

                    for iloop_MOD=1:length(set_MOD)
                        MODCF_type = set_MOD{iloop_MOD}

                        for iloop_abs=1:length(set_real_abs)
                            error_type= set_real_abs{iloop_abs};

                            icount_DRIVER = icount_DRIVER+1;


                            ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
                            modis_data_plot = modis_data_plots{iloop_plot};
                            
                            switch modis_data_plot
                                case 'Number of droplets cell values time mean - specific days'
                                    iset_min_clim=0;
                                    clim_min=0;
                                    iset_max_clim=0;
                                    %        clim_max=225;
                                    clim_max=150;
                                otherwise
                                    iset_min_clim=1;
                                    clim_min=-20;
                                    iset_max_clim=1;
                                    %        clim_max=225;
                                    clim_max=20;

                            end

                            icontour = icontour_set;
                                                        
                            plot_global_maps
                            Psave{icount_DRIVER} = P;
                            Npoints_save_loop{icount_DRIVER} = Npoints;
                            threshCF_save_loop{icount_DRIVER} = thresh_CF;
                            
                            if isave_vals==1 & exist('mean_bias_MODIS')
                                mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
                                RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;

                                fprintf(fid_mtab,'%f ',mean_bias_MODIS);
                                fprintf(fid_Rtab,'%f ',RMSE_MODIS);
                            end
                          
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
%do a second plot of the fraction of days that were left after screening  -- PROB don't need this.                       
                             ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
                            modis_data_plot = fraction_plotcase;
                            if fraction_days_plot ==1
                                plot_global_maps
                                if isave_plot==1
                                    saveas_ps_fig_emf(gcf,[savename],'',0,1);
                                end
                            end
                          



                        end

                    end

                end

            end

        end
        
%              save(save_file,'modis_data_plots','set_screening','set_MOD','set_real_abs','Psave','mean_bias_DRIVER','RMSE_DRIVER','screening_eval_str','screen_edits')
           
        if isave_vals==1
            fprintf(fid_mtab,'\n');
            fprintf(fid_Rtab,'\n');
        end

     
        
        if isave_vals==1
            fclose(fid_mtab);
            fclose(fid_Rtab);
        end
        
        
%Now plot the difference

colormap_choose=jet; %default

            P = 100*(Psave{1}-Psave{2})./Psave{2};
            Npoints = min(Npoints_save_loop{1},Npoints_save_loop{2});
                
            MODIS_varname2_plot = ['Percentage diff CF>' num2str(threshCF_save_loop{1}(1)) '.AND.LT.' num2str(threshCF_save_loop{1}(2)) ' screening vs all CF for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str];

                units_str_plot='%';
                title_info = '';
                icontour=0;
                inew_cticks=0;


                iset_min_clim=1;
                clim_min=-25;
                iset_max_clim=1;
                clim_max=25;
                
                ifilter_ndays=0;
                thresh_ndays=15;
        
                
                ioverride_plotglobal_thresh=1;                
                ioverride_time_selection=1;
                ioverride_plotglobal_loc=1;
                mod_data_type='timeseries3 lambert';
                modis_data_plot = 'Generic plot specified outside of script';

                plot_global_maps
                if isave_plot==1
                    saveas_ps_fig_emf(gcf,[savename],'',0,1);
                end
        
                          

                            
        
 case 'MODIS vs AMSRE LWP bias with Ndays and mean CF contours'    
        %runs through different screenings, shows the mean CF contours in
        %order to assess impact of CF screening on amount of data excluded.
       %Also plots the fraction of possible days used
       
       fraction_days_plot = 0; %flag for wheterh to do the fraction of days remaining after screening plot
       
        modis_data_plots={'MODIS LWP minus AMSRE'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';
%        modis_data_plots={'Number of droplets cell values time mean - specific days'}; fraction_plotcase='Fraction of days remaining after screening';
%         modis_data_plots={'LWP AMSRE time3'}; fraction_plotcase='Fraction of days remaining after screening for AMSRE bias comp';

%        zero_CF_out_of_range=0;
        time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';
        
        set_screening = {'none','NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'};
        % --- NOTE the below uses CTH_max, not CTH ---
        set_screening = {'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff'};


%set a default screening - then can change individual ones in the loop
screening_eval_str=[ ...
        'thresh_ndays=15;' ...
        'thresh_SZA=[0 90];' ...
        'thresh_CF=[0.8 1.00001];' ...
        'thresh_NP=50;' ...
        'thresh_sensZA=[0 90];'     ...
        'thresh_CTT = [273 273+100];'  ...
        'thresh_relAZ = [0 180];'   ...      
        'thresh_CTH = [-0.01 3.2];' ...
        'thresh_stdW = [0 1e9];' ...
        'thresh_sigCTT = [0 1e9];' ...
        'thresh_reff = [0 14];' ...
        ];

    %just ones that alter it from the above baseline
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];'};
screen_edits = {'thresh_CF=[0.8 1.00001]','thresh_CF=[0.5 1.00001]','thresh_CF=[-0.01 1.00001]','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9]; thresh_reff=[0 30];','thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];','thresh_CF=[0.8 1.00001]; thresh_CTH=[-0.01 1e9];','thresh_CF=[-0.01 1.00001]; thresh_CTH=[-0.01 1e9];'};
screen_edits = {'thresh_reff=[0 30];'};
screen_edits = {''};
%screen_edits = {'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];'};
screen_edits = {'','thresh_CF=[-0.01 1.00001];', 'thresh_CF=[-0.01 1.00001]; thresh_reff=[0 30];',...
    'thresh_CF=[0.1 0.4];' };


        
       
        
        
        set_MOD = {'MOD06 liquid','MOD35'};
%         set_MOD = {'MOD35'};
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
        set_real_abs = {'absolute','percentage'};
        set_real_abs = {'percentage'};   
%        set_real_abs = {'absolute'};
        

        ifilter_ndays=0;  %filtering for number of data points
        ifilter_clhcalipso=0; %whether to do the filtering using CALIPSO high cloud or not
        
        %what to plot for the black contours.
        cont_dat_choose = 'calipso mid+highCF';
        cont_col_str = 'k';
        
        %what to use for filtering
        filt_dat_choose = 'calipso mid+highCF';
        thresh_clh = [-0.01 0.3];close all
        
        iuse_saved_mask=0; %This uses a saved mask to make NaNs - might override the ifilter_clhcalipso
        saved_mask_name = 'imask_ndays15_cal_midhigh';
        
        clear RMSE_DRIVER mean_bias_DRIVER Psave
        icount_DRIVER=0;             
        
        for iloop_plot=1:length(modis_data_plots)  %length=1
            modis_data_plot = modis_data_plots{iloop_plot};

            for iscreen_set=1:length(set_screening) %length=1
                screen_type = set_screening{iscreen_set};

                for iscreen_thresh=1:length(screen_edits) %length is variable
                    eval(screening_eval_str); %eval the basline screening
                    eval(screen_edits{iscreen_thresh}); %and now the edits to the baseline

                    for iloop_MOD=1:length(set_MOD) %MOD35 and MOD06
                        MODCF_type = set_MOD{iloop_MOD}

                        for iloop_abs=1:length(set_real_abs)  %generally just set for 'percentage'
                            error_type= set_real_abs{iloop_abs};

                            icount_DRIVER = icount_DRIVER+1;


                            ioverride_plotglobal_thresh=1;
                            iocean_only=0;
                            iadd_correct_for_clear_sky_bias=0; %Subtract the clear-sky bias for AMSRE (done in plot_global_maps)
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
                            modis_data_plot = modis_data_plots{iloop_plot};
                            
                            switch modis_data_plot
                                case 'LWP AMSRE time3'
                                    iset_min_clim=1;
                                    clim_min=0;
                                    iset_max_clim=1;
                                    %        clim_max=225;
                                    clim_max=150;
                                otherwise
                                    iset_min_clim=1;
                                    clim_min=-20;
                                    iset_max_clim=1;
                                    %        clim_max=225;
                                    clim_max=20;

                            end

                            
                                                        
                            plot_global_maps
                            Psave{icount_DRIVER} = P;
%                            mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
%                            RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;
   %Better to use these as they have the filtering for ndays taken into
   %account (and the domain restriction)
                            mean_bias_DRIVER(icount_DRIVER)=Pmean
                            RMSE_DRIVER(icount_DRIVER)=P_RMSE                            
                            
                            if isave_vals==1
                                fprintf(fid_mtab,'%f ',Pmean);
                                fprintf(fid_Rtab,'%f ',P_RMSE);
                            end
                          
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
%do a second plot of the fraction of days that were left after screening                            
                             ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='timeseries3 lambert';
                            modis_data_plot = fraction_plotcase;
                            if fraction_days_plot ==1
                                plot_global_maps
                                if isave_plot==1
                                    saveas_ps_fig_emf(gcf,[savename],'',0,1);
                                end
                            end
                          



                        end

                    end

                end

            end

        end
           
        if isave_vals==1
            fprintf(fid_mtab,'\n');
            fprintf(fid_Rtab,'\n');
        end

     
        
        if isave_vals==1
            fclose(fid_mtab);
            fclose(fid_Rtab);
        end
        
        
%       save(save_file,'modis_data_plots','set_screening','set_MOD','set_real_abs','Psave','mean_bias_DRIVER','RMSE_DRIVER','screening_eval_str','screen_edits')
        
                    eval(screening_eval_str); %eval the basline screening
                    eval(screen_edits{iscreen_thresh}); %and now the edits to the baseline
        
        
    case 'MODIS vs AMSRE LWP bias'    
        %runs through different screenings and also does MOD35 and MOD06 vs
        %AMSRE
        modis_data_plots={'MODIS LWP minus AMSRE'}; 
%        zero_CF_out_of_range=0;
        time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';
        
        set_screening = {'none','NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW'};
        set_screening = {'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW'};
        thresh_ndays=100; %threshold no. days 
        thresh_SZA=[0 90];
        thresh_CF=[0.8 1.00001];
        thresh_NP=50;
        thresh_sensZA=[0 90];    
        thresh_CTT = [273 273+100]; 
        thresh_relAZ = [0 180];         
        thresh_CTH = [-0.01 3.2];
        thresh_stdW = [0 1e9];
        thresh_sigCTT = [0 1e9];
        
        
        
        
        
        ifilter_ndays=1;
        
        iset_min_clim=1;
        clim_min=20;
        iset_max_clim=1;
%        clim_max=225;
        clim_max=150;
        
        
        set_MOD = {'MOD06 liquid','MOD35'};
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
        set_real_abs = {'absolute','percentage'};
        
        clear RMSE_DRIVER mean_bias_DRIVER Psave
        icount_DRIVER=0;             
        
        for iloop_plot=1:length(modis_data_plots)
            modis_data_plot = modis_data_plots{iloop_plot};

            for iscreen_set=1:length(set_screening)
                screen_type = set_screening{iscreen_set};              
            
                  for iloop_MOD=1:length(set_MOD)
                      MODCF_type = set_MOD{iloop_MOD}
                      
                      for iloop_abs=1:length(set_real_abs)
                          error_type= set_real_abs{iloop_abs};
                      
                      icount_DRIVER = icount_DRIVER+1;
                      
            
                      ioverride_plotglobal_thresh=1;
%                      iocean_only=1;
                      ioverride_time_selection=1;
                      ioverride_plotglobal_loc=1;                      
                      mod_data_type='timeseries3 lambert';
                      plot_global_maps
                      Psave{icount_DRIVER} = P;
                      mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;                      
                      RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;
                      
                   if isave_vals==1   
                      fprintf(fid_mtab,'%f ',mean_bias_MODIS);
                      fprintf(fid_Rtab,'%f ',RMSE_MODIS);                      
                   end
                   
                      if min(thresh_clh==[-0.01 1.01])==1
                          %save if we plotted the one with no data clipping
                          %using CALIPSO CF
                          
%                          saveas_ps_fig_emf(gcf,[savename],'',0,1);
                      end
                      
                       if isave_plot==1
                          saveas_ps_fig_emf(gcf,[savename],'',0,1);
                      end

                  end

            end
        
            end
        
        end
           
        if isave_vals==1
            fprintf(fid_mtab,'\n');
            fprintf(fid_Rtab,'\n');
        end

     
        
        if isave_vals==1
            fclose(fid_mtab);
            fclose(fid_Rtab);
        end
        
%        save(save_file,'modis_data_plots','set_screening','set_MOD','set_real_abs','Psave','mean_bias_DRIVER','RMSE_DRIVER')
        

    case 'MODIS vs AMSRE LWP'    
        modis_data_plots={'LWP AMSRE time3','LWP cell values time mean (grid-box mean using MOD35) - specific days'}; 
%        zero_CF_out_of_range=0;
        time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';
        
        set_screening = {'none','NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW'};
%        set_screening = {'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW'};
        thresh_ndays=15; %threshold no. days 
        thresh_SZA=[0 90];
        thresh_CF=[0.8 1.00001];
        thresh_NP=50;
        thresh_sensZA=[0 90];    
        thresh_CTT = [273 273+100]; 
        thresh_relAZ = [0 180];         
        thresh_CTH = [-0.01 3.2];
        thresh_stdW = [0 1e9];
        thresh_sigCTT = [0 1e9];
        
        ifilter_ndays=1;
        
        iset_min_clim=1;
        clim_min=20;
        iset_max_clim=1;
%        clim_max=225;
        clim_max=150;
        
        
%        set_MOD = {'MOD06 liquid','MOD35'};
%        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};
        
        clear RMSE_DRIVER mean_bias_DRIVER
        icount_DRIVER=0;             
        
        for iloop_plot=1:length(modis_data_plots)
            modis_data_plot = modis_data_plots{iloop_plot};

            for iscreen_set=1:length(set_screening)
                screen_type = set_screening{iscreen_set};              
            
                  
                      
                      icount_DRIVER = icount_DRIVER+1;
                      
            
                      ioverride_plotglobal_thresh=1;
%                      iocean_only=1;
                      ioverride_time_selection=1;
                      ioverride_plotglobal_loc=1;                      
                      mod_data_type='timeseries3 lambert';
                      plot_global_maps
                      mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;                      
                      RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;
                      
                   if isave_vals==1   
                      fprintf(fid_mtab,'%f ',mean_bias_MODIS);
                      fprintf(fid_Rtab,'%f ',RMSE_MODIS);                      
                   end
                   
                      if min(thresh_clh==[-0.01 1.01])==1
                          %save if we plotted the one with no data clipping
                          %using CALIPSO CF
                          
%                          saveas_ps_fig_emf(gcf,[savename],'',0,1);
                      end
                      
                      if isave_plot==1
                          saveas_ps_fig_emf(gcf,[savename],'',0,1);
                      end
                      
                  end

        end
           
        if isave_vals==1
            fprintf(fid_mtab,'\n');
            fprintf(fid_Rtab,'\n');
        end

     
        
        if isave_vals==1
            fclose(fid_mtab);
            fclose(fid_Rtab);
        end

        
        
case 'MODIS vs CALIPSO bias'
    
        icalc_saved_modis_monthly_CFs=0; %Flag whether to calculate the monthly MODIS CFs for the different screenings, or whether
            %to load in a saved .mat file
    
        modis_data_plot='MODIS Cloud Fraction minus CALIPSO CF'; mod_data_type='timeseries3 lambert';
%        zero_CF_out_of_range=0;
        time_mean_str ='ALL'; %This chooses the times required
        cont_dat_choose = 'calipso mid+highCF';    
    


         modis_calipso_cf_case = 'calc and save MODIS monthly data'; %Compute and save the MODIS monthly CF data for the different screenings
        modis_calipso_cf_case = 'mean time bias'; %the old case where do the time mean to get one value for each location
%        modis_calipso_cf_case = 'RMSE from monthly data'; %Doing RMSE on monthly data as a map for each location
%        modis_calipso_cf_case = 'Correlation coeff from monthly data'; %Doing RMSE on monthly data as a map for each location        
        
%        modis_save_monthly_file = '~/modis_saved_monthly_2006_2009_data.mat';
%        modis_save_monthly_file = '~/modis_saved_monthly_2007_2010_data.mat';        
        modis_save_monthly_file = '~/modis_saved_monthly_2007_2010_data_031616.mat';              
%        modis_save_monthly_file = '~/modis_saved_monthly_2007_2010_data_2deg_avs_031716.mat';                      
        
        iclose_figs_cal_MODIS_cf=0; %Gets set to one for the monthly calculation
        
        
%% Set the screening to be used in modisL3_screening_timeseries3       
        % Screening used in original draft of the paper was 'CTH only with
        % no zeroCF screening' - showed tests with and without the CF set
        % to zero instead of NaN when out of range (set_screening indices 2 and 3
        % below).
        
        %Order of the loops are set_thresh_CTH, set_calCF, set_MOD and set_screening from
        %outer to inner.
        
%        set_izero = {0,0,1,0,0,0,0,0}; %Whether to set zero_CF_out_of_range (setting to zero instead of NaN test)
            %these flags correspond to set_screening below, so make the
            %same length

%        set_screening = {'none','CTH only with no zeroCF screening','CTH only with no zeroCF screening','CTP only with no zeroCF screening',...
%            'CTP and CTH only with no zeroCF screening','Min CTP only with no zeroCF screening','max CTH only with no zeroCF screening',...
%            'min_CTP and max_CTH hybrid only with no zeroCF screening'};

 set_izero = {0,0,0,0,0}; %Whether to set zero_CF_out_of_range (setting to zero instead of NaN)
            %these flags correspond to set_screening below, so make the
            %same length
% Cutting down these screenings        
        set_screening = {'none','CTH only with no zeroCF screening',...
            'min_CTP and max_CTH hybrid only with no zeroCF screening','CTH_Ryan only with no zeroCF screening','CTH_adjust only with no zeroCF screening'};
                
        
        set_thresh_CTH = {[-0.01 3.2],[-0.01 2]}; %N.B. - this makes no difference to the hybrid screening because the height at which
          % the screening switches to CTP is set in filtering_dat_get (to
          % 2km at the moment).
        
        thresh_CTH = [-0.01 3.2];
%        thresh_CTH = [-0.01 2];  %Reducing to 3.2km doesn't make too much
%                                    difference either...
%        thresh_CTH = [-0.01 6]; %Increasign made it worse...
%        thresh_CTH = [-0.01 1.0];  %Reducing to 1.5 or 1km makes better...
                                  
        thresh_CTP = [680 1200]; %Cloud top pressure screening - relying on the CO2 slicing being mostly ok for CTP<700 hPa or so
        set_MOD = {'MOD06 liquid','MOD35'};

        set_calCF = {[-0.01 1.01],[-0.01 0.3],[0.3 1.01]};  %Loops over all these for each screening and MOD35/MOD06 setting
        

        
        switch modis_calipso_cf_case                
            case 'calc and save MODIS monthly data'
                imonths_driver = [1:12];
                iyears_driver = [2007:2010];
                ioverride_years_time_screen=1;
                modisyear_timeseries3 = modisyear_timeseries3_bk; %reset this since it may have got overwritten
                siz=size(cllcalipso_monthly);
                siz_modis_cf_save = [length(iyears_driver) siz(2) siz(3) length(set_thresh_CTH) length(set_calCF) length(set_MOD) length(set_screening)];
                %Empty array for saving the monthly screened MODIS data interpolated
                %onto the CALIPSO grid
                modis_cf_grid_save = NaN*ones( siz_modis_cf_save );
                iclose_figs_cal_MODIS_cf=1;
                
            case {'mean time bias','RMSE from monthly data','Correlation coeff from monthly data'}
                imonths_driver = [0];
                iyears_driver = [0];
                switch modis_calipso_cf_case
                    case 'RMSE from monthly data'
                        load(modis_save_monthly_file);
                end
        end
        
  
        
        clear RMSE_DRIVER mean_bias_DRIVER corr_coeff_DRIVER
        icount_DRIVER=0;
        
   for iyear_modis_vs_calipso=1:length(iyears_driver)   
       
       years_required_for_mean = iyears_driver(iyear_modis_vs_calipso); % for time_inds_modisL3_timeseries3
        
     for imonth_modis_vs_calipso=1:length(imonths_driver)  
         
         switch modis_calipso_cf_case
             case 'calc and save MODIS monthly data'
                %get the day of year for the month required
                [days_required_for_mean,time_mean_str] = days_of_month(imonth_modis_vs_calipso);
                time_index_calipso = (iyear_modis_vs_calipso-1)*12 + imonth_modis_vs_calipso;                 
             otherwise
                 time_index_calipso=[1:48]; %all times
         end
         
      for iCTH_set=1:length(set_thresh_CTH)
          thresh_CTH = set_thresh_CTH{iCTH_set};
        
        for ical_CF = 1:length(set_calCF) 
            ifilter_clhcalipso=1;
            thresh_clh =  set_calCF{ical_CF};

        
        for iMOD_set=1:length(set_MOD) 
            MODCF_type = set_MOD{iMOD_set};

            for iscreen_set=1:length(set_screening)
                screen_type = set_screening{iscreen_set};
                zero_CF_out_of_range = set_izero{iscreen_set};

                switch screen_type
                    case 'none'              
                        cal_CF_range='low+mid+high';
%                    case 'CTH only with no zeroCF screening'
                    otherwise
                        cal_CF_range='low';
                end

            

                  
                      
                      icount_DRIVER = icount_DRIVER+1;
                      
            
                      ioverride_plotglobal_thresh=1;
                      iocean_only=1;
                      ioverride_time_selection=1;
                      ioverride_plotglobal_loc=1;
                      plot_global_maps
                      
                      if iclose_figs_cal_MODIS_cf==1
                        close(gcf);
                      end
                      
                      mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;                      
                      RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;
                      corr_coeff_DRIVER(icount_DRIVER)=RHO;
                      
                      score_corr_DRIVER(icount_DRIVER)=score_from_one;
                      
                   if isave_vals==1   
                      fprintf(fid_mtab,'%f ',mean_bias_MODIS);
                      fprintf(fid_Rtab,'%f ',RMSE_MODIS);                      
                   end
                   
                      if min(thresh_clh==[-0.01 1.01])==1
                          %save if we plotted the one with no data clipping
                          %using CALIPSO CF
                          
%                          saveas_ps_fig_emf(gcf,[savename],'',0,1);
                      end

                  end

        end
           
        if isave_vals==1
            fprintf(fid_mtab,'\n');
            fprintf(fid_Rtab,'\n');
        end

        end
        
     end
     
     end
end
        
        if isave_vals==1
            fclose(fid_mtab);
            fclose(fid_Rtab);
        end
        
%       save(modis_save_monthly_file,'modis_cf_grid_save','set_thresh_CTH','set_calCF','set_MOD','set_screening','set_izero');
        
      
        %Choose some screenings to output
         %Order of the loops are set_calCF, set_MOD and set_screening from
        %outer to inner.
        

        L_cal = length(set_calCF);
        L_mod = length(set_MOD);
        L_screen = length(set_screening);
        
        iCTHs=[1 2];
        icals=[3]; % 1=all regions, 2=low CF region only, 3=high CF region only
        iscreens=[1 2 3];

        
        fprintf(1,'\n');
        fprintf(1,'  i)     Overall RMSE    Corr coeff    Score vs r=1 \n');
        for iCTH=iCTHs
            for ical=icals
                for imod=[1:2]
                    for iscreen=iscreens
                        i = (iCTH-1)*L_cal*L_mod*L_screen + (ical-1)*L_mod*L_screen + (imod-1)*L_screen + iscreen;
                        fprintf(1,'\n%d    %.3f    %.3f    %.3f',i,RMSE_DRIVER(i),corr_coeff_DRIVER(i),score_corr_DRIVER(i));
                    end
                end
            end
        end
        fprintf(1,'\n');
        
        
    case 'GCM CFs'

%DAYTIME
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
time_mean_str ='ANNUAL';
modis_data_plot='Max Low Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
times_required = [12:15];  %Aqua daytime
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);

ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
time_mean_str ='ANNUAL';
modis_data_plot='Max Mid Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
times_required = [12:15];  %Aqua daytime
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);
 
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
time_mean_str ='ANNUAL'; 
modis_data_plot='Max High Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
times_required = [12:15];  %Aqua daytime
plot_global_maps                 
saveas_ps_fig_emf(gcf,[savename],'',0,1);

%NIGHTTIME
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
time_mean_str ='ANNUAL';
modis_data_plot='Max Low Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
times_required = [0:3];  %Aqua nighttime
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);

ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
time_mean_str ='ANNUAL';
modis_data_plot='Max Mid Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
times_required = [0:3];  %Aqua nighttime
plot_global_maps
saveas_ps_fig_emf(gcf,[savename],'',0,1);
 
ioverride_plotglobal_thresh=1;
ioverride_time_selection=1;
time_mean_str ='ANNUAL'; 
modis_data_plot='Max High Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
times_required = [0:3];  %Aqua nighttime
plot_global_maps                 
saveas_ps_fig_emf(gcf,[savename],'',0,1);

end


clear ioverride_DRIVER
catch error_DRIVER
    clear ioverride_DRIVER
    rethrow(error_DRIVER)
end
