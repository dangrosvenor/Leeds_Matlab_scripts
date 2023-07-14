       % Creates maps of mean difference in LWP between model and AMSRE
       % (based on the time-means for each location).
       % Taking the approach of running plot_global once for MODIS and then once
       % for the model since they use different screenings. Then plotting
       % the difference.
       
       iadd_correct_for_clear_sky_bias = 1; %Whether to correct AMSRE for the daytime clear-sky  bias
             %as calculated using MODIS CF
       
       fraction_days_plot = 0; %flag for wheterh to do the fraction of days remaining after screening plot

       LWP_thresh_choose = 50; %ignore any data with LWP less than this (g/m2).       
       LWP_thresh_choose = -1e9; %ignore any data with LWP less than this (g/m2).
       
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
        set_screening = {'none'};   %Can't use MODIS screenings here because want to compare at night time too
                                    %So, prob don't need the MODIS data
                                    %loaded??
                                    
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
        

        
        
        amsre_daynight_labels = {'daytime','nighttime'}; %loop thorugh day and night plots

        
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

                    for idaynight=1:length(amsre_daynight_labels)
                        amsre_daynight_label= amsre_daynight_labels{idaynight};

                        for iloop_abs=1:length(set_real_abs)
                            error_type= set_real_abs{iloop_abs};
                            icount_DRIVER = icount_DRIVER+1;
                            
                    % -- run AMSRE LWP first        
%                            modis_data_plot ='Number of droplets cell values time mean - specific days';
                            modis_data_plot ='LWP AMSRE conditional sampling'; 
                            gcm_str_select = 'AMSRE';
                            




                            ioverride_plotglobal_thresh=1;
                            %                      iocean_only=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='AMSRE';
                            LWP_thresh = LWP_thresh_choose; %ignore any data with LWP less than this (g/m2).                            
%                            modis_data_plot = modis_data_plots{iloop_plot};
                            

                                    iset_min_clim=0;
                                    clim_min=0;
                                    iset_max_clim=0;
                                    %        clim_max=225;
                                    clim_max=150;
                                    
                                    
%                                    amsre_daynight_label= 'daytime';
%                                    amsre_daynight_label= 'nighttime';
  %set above in loop using amsre_daynight_label= amsre_daynight_labels{idaynight};

                                    
%                            times_required = [0:24]; %need to specify all hours requried (in 1 hour increments)
%                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                            years_required_for_mean = [2006:2010]; time_mean_str = 'ANNUAL 2006-2010';
                            months_required_for_mean = [1:12];
                               %Need to change time_mean_str to somethign
                               %different to 'ALL' or 'ANNUAL' as will
                               %think it's MODIS data
                            ioverride_years_time_screen=1; %Also need to set this
                                                        
                            plot_global_maps
                            Psave_DRIVER{icount_DRIVER} = P_save;
                            Psave_MODIS = P_save;                            
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
                    % -- now do the model LWP    
                            modis_data_plot ='LWP GCM (grid-box mean, no CF screening) conditional sampling'; mod_data_type='GCM';   gcm_time_of_day_select=2;                
                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
%                            years_required_for_mean = [2006:2010]; days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL 2006-2010';

                    switch amsre_daynight_label
                        case 'daytime' 
                            times_required = [12:15];  %Aqua daytime
                        case 'nighttime'
                            times_required = [0:3]; %Aqua nighttime
                    end

                            
                            gcm_str_select = gcm_str_last_loaded;
                            
                            %tells 'Max Nd in low cloud without cloud screening' what screening to do.
                            screen_using_cosp = 're CF range'; 
                            screen_using_cosp = 're';
                            screen_using_cosp = 'none';                            

                            LWP_thresh = LWP_thresh_choose; %ignore any data with LWP less than this (g/m2).                            

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
                            
                                    icount_DRIVER = icount_DRIVER+1;
                            plot_global_maps
                            %interpolate the model field to the AMSRE grid
                            %Uses P_save because this is the P field before
                            %any scaling for the colobar etc.
%                            P_int =  griddata(Plon2D,Plat2D,P_save,MLON,MLAT');  %for MODIS
                            P_int =  griddata(Plon2D,Plat2D,P_save,gcm_Plon2D_AMSRE,gcm_Plat2D_AMSRE); %for AMSRE
                            Psave_DRIVER{icount_DRIVER} = P_int;
                            Psave_model = P_int;
                            Npoints_save_model = Npoints;
%                            mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
%                            RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;

                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end                                                                                    

                            
                            
                            
%  ---               Now plot the difference

%                            colormap_choose=jet; %default
                            lb_map=lbmap(16,'brownblue'); %nice colormap for colorblind people
        
                             i_reverse_cmap=1;
                             if i_reverse_cmap==1
                                 cmap=flipud(lb_map);
                             else
                                 cmap=lb_map;
                             end
                             colormap_choose=cmap;


                             P = 100*(Psave_model-Psave_MODIS)./Psave_MODIS;
                             Npoints = min(Npoints_save_model,Npoints_save_MODIS);

%                            MODIS_varname2_plot = ['Percentage diff in N_d for model (' gcm_str ' using ' modis_data_plot cosp_screen_str ') vs MODIS for ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str_MODIS];
                            MODIS_varname2_plot = ['Percentage diff bias for model (' gcm_str ' using ' modis_data_plot ') vs AMSRE for ' time_mean_str ' ' thresh_str_MODIS];                            

                            units_str_plot='%';
                            title_info = '';
                            
                            %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=0;
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
                            mod_data_type='AMSRE';
                            gcm_str_select = 'AMSRE';
                            modis_data_plot = 'Generic plot specified outside of script';

                                                        
%                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                            years_required_for_mean = [2006:2010]; time_mean_str = 'ANNUAL 2006-2010';
                            months_required_for_mean = [1:12];
                               %Need to change time_mean_str to somethign
                               %different to 'ALL' or 'ANNUAL' as will
                               %think it's MODIS data
                            ioverride_years_time_screen=1; %Also need to set this
                            
                            icount_DRIVER = icount_DRIVER+1;
                            plot_global_maps
                            Psave_DRIVER{icount_DRIVER} = P_save;
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
                            if isave_vals==1
                                fprintf(fid_mtab,'%f ',Pmean);
                                fprintf(fid_Rtab,'%f ',P_RMSE);
                            end



                        end

                    end

                end

            end

        end
        
        
        %  ---               Now plot the diurnal range bias

%                            colormap_choose=jet; %default
                            lb_map=lbmap(16,'brownblue'); %nice colormap for colorblind people
        
                             i_reverse_cmap=1;
                             if i_reverse_cmap==1
                                 cmap=flipud(lb_map);
                             else
                                 cmap=lb_map;
                             end
                             colormap_choose=cmap;
                             
                            %In the loop have saved the AMSRE value, the
                            %model value and the difference (bias in
                            %model). Does this for both day and night. So
                            %there are 6 saved P fields.
                            diurnal_range_AMSRE = Psave_DRIVER{4} - Psave_DRIVER{1};
                            diurnal_range_model = Psave_DRIVER{5} - Psave_DRIVER{2};
                            P = 100*(diurnal_range_model-diurnal_range_AMSRE)./diurnal_range_AMSRE;
                            Npoints = min(Npoints_save_model,Npoints_save_MODIS);

%                            MODIS_varname2_plot = ['Percentage diff in N_d for model (' gcm_str ' using ' modis_data_plot cosp_screen_str ') vs MODIS for ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str_MODIS];
                            MODIS_varname2_plot = ['Percentage diff in diurnal LWP range for model (' gcm_str ' using ' modis_data_plot ') vs AMSRE for ' time_mean_str ' ' thresh_str_MODIS];                            

                            units_str_plot='%';
                            title_info = '';
                            
                            %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=0;
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
                            mod_data_type='AMSRE';
                            gcm_str_select = 'AMSRE';
                            modis_data_plot = 'Generic plot specified outside of script';

                                                        
%                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                            years_required_for_mean = [2006:2010]; time_mean_str = 'ANNUAL 2006-2010';
                            months_required_for_mean = [1:12];
                               %Need to change time_mean_str to somethign
                               %different to 'ALL' or 'ANNUAL' as will
                               %think it's MODIS data
                            ioverride_years_time_screen=1; %Also need to set this
                            
                            plot_global_maps
                            Psave_DRIVER{icount_DRIVER} = P_save;
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            
                            if isave_vals==1
                                fprintf(fid_mtab,'%f ',Pmean);
                                fprintf(fid_Rtab,'%f ',P_RMSE);
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
        
   

%  ---               Also plot the diurnal range (rather than the bias as
%  done above)

% ------------------  Model diurnal range -------------------------

%                            colormap_choose=jet; %default
                            lb_map=lbmap(16,'brownblue'); %nice colormap for colorblind people
        
                             i_reverse_cmap=1;
                             if i_reverse_cmap==1
                                 cmap=flipud(lb_map);
                             else
                                 cmap=lb_map;
                             end
                             colormap_choose=cmap;
                             
                            %In the loop have saved the AMSRE value, the
                            %model value and the difference (bias in
                            %model). Does this for both day and night. So
                            %there are 6 saved P fields.
%                            diurnal_range_AMSRE = Psave_DRIVER{4} - Psave_DRIVER{1};
%                            diurnal_range_model = Psave_DRIVER{5} - Psave_DRIVER{2};
                            P = diurnal_range_model;
                            Npoints = min(Npoints_save_model);

%                            MODIS_varname2_plot = ['Percentage diff in N_d for model (' gcm_str ' using ' modis_data_plot cosp_screen_str ') vs MODIS for ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str_MODIS];
                            MODIS_varname2_plot = ['Diurnal LWP range for model (' gcm_str_last_loaded ' using ' modis_data_plot ') for ' time_mean_str];                            

                            units_str_plot='g m^{-2}';
                            title_info = '';
                            
                            %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=0;
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
                            clim_min=0;
                            iset_max_clim=1;
                            clim_max=100;

                            ifilter_ndays=0;
                            thresh_ndays=15;


                            ioverride_plotglobal_thresh=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='AMSRE';
                            gcm_str_select = 'AMSRE';
                            modis_data_plot = 'Generic plot specified outside of script';

                                                        
%                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                            years_required_for_mean = [2006:2010]; time_mean_str = 'ANNUAL 2006-2010';
                            months_required_for_mean = [1:12];
                               %Need to change time_mean_str to somethign
                               %different to 'ALL' or 'ANNUAL' as will
                               %think it's MODIS data
                            ioverride_years_time_screen=1; %Also need to set this
                            
                            plot_global_maps
                            Psave_DRIVER{icount_DRIVER} = P_save;
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                            

%  ---       AMSRE dirunal range --------------------------

%                            colormap_choose=jet; %default
                            lb_map=lbmap(16,'brownblue'); %nice colormap for colorblind people
        
                             i_reverse_cmap=1;
                             if i_reverse_cmap==1
                                 cmap=flipud(lb_map);
                             else
                                 cmap=lb_map;
                             end
                             colormap_choose=cmap;
                             
                            %In the loop have saved the AMSRE value, the
                            %model value and the difference (bias in
                            %model). Does this for both day and night. So
                            %there are 6 saved P fields.
%                            diurnal_range_AMSRE = Psave_DRIVER{4} - Psave_DRIVER{1};
%                            diurnal_range_model = Psave_DRIVER{5} - Psave_DRIVER{2};
                            P = diurnal_range_AMSRE;
                            Npoints = Npoints_save_MODIS;

%                            MODIS_varname2_plot = ['Percentage diff in N_d for model (' gcm_str ' using ' modis_data_plot cosp_screen_str ') vs MODIS for ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str_MODIS];
                            MODIS_varname2_plot = ['Diurnal LWP range for AMSRE for ' time_mean_str ' ' thresh_str_MODIS];                            

                            units_str_plot='g m^{-2}';
                            title_info = '';
                            
                            %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=0;
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
                            clim_min=0;
                            iset_max_clim=1;
                            clim_max=100;

                            ifilter_ndays=0;
                            thresh_ndays=15;


                            ioverride_plotglobal_thresh=1;
                            ioverride_time_selection=1;
                            ioverride_plotglobal_loc=1;
                            mod_data_type='AMSRE';
                            gcm_str_select = 'AMSRE';
                            modis_data_plot = 'Generic plot specified outside of script';

                                                        
%                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                            years_required_for_mean = [2006:2010]; time_mean_str = 'ANNUAL 2006-2010';
                            months_required_for_mean = [1:12];
                               %Need to change time_mean_str to somethign
                               %different to 'ALL' or 'ANNUAL' as will
                               %think it's MODIS data
                            ioverride_years_time_screen=1; %Also need to set this
                            
                            plot_global_maps
                            Psave_DRIVER{icount_DRIVER} = P_save;
                            if isave_plot==1
                                saveas_ps_fig_emf(gcf,[savename],'',0,1);
                            end
                                                        