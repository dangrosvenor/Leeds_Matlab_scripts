try
    
    icontour_DRIVER=0;
    igeneric_plot=0;
    
    min_LWP = 1.0; %g/m2
    min_CF = 0.05; %min CF allowed to avoid divide by zero for in-cloud LWP calc

    isub=0; clear hs %for subplots
    
    if ~exist('ioverride_LWP_Sc_choices') | ioverride_LWP_Sc_choices==0
        
        
                        
        % Which variable to plot? :-
        plot_var_Sc = 'LWP'; %Plot the LWP and bias
        plot_var_Sc = 'LWPic'; %Plot the LWP divided by the CALIPSO low CF, and bias
        %plot_var_Sc = 'Ndat'; % Plot the number of days that pass the filter criteria, or the LWP
        %plot_var_Sc = 'Ratio of LWP to TLWP'; %
        %plot_var_Sc = 'COSP MODIS low CF estimate test';
        %plot_var_Sc = 'Nd';
        %plot_var_Sc = 'Nd COSP';
      
        
        bias_type = 'Absolute';
        bias_type = 'Percentage'; %overrides to absolute if plotting ratio of LWP to TLWP
                
        ibias_contour=0;
        
        %----------------------------------------------------------------
        % Switches that affect both MODEL and SATELLITE
        %----------------------------------------------------------------
        
        %Set which AMSRE data to compare to. If want to match to the MODIS data to do the MODIS-based screening (low clouds only, etc.)
        %then need to set to 'day'. Then will also match the MODEL to
        %aorund 13 LT. Night will match the model to 01 LT. Average will
        %use an average of both.            
        lwp_day_night = 'day';
        lwp_day_night = 'night';        
        lwp_day_night = 'average';
        
        iadd_RWP=1; %whether to include RWP as well as LWP (model and sat)
        
        irestrict_using_LWP_fraction=0; %whether to restrict to boxes with a high LWP to TLWP ratio
        min_LWP_ratio=0.99; 
        
        iNtot_ratio_CF80_to_0 = 1; %flag that sets method for calculating the ratio of Sc days to the total. Setting to one
        %uses the ratio of CF>80 to CF>=0 clouds. Clouds are still restricted
        %to low alt clouds by using mean CTH<=3.2km for MODIS and using only
        %gridboxes with no mid or high level cloud in the model (using COSP
        %MODIS low, mid and high clouds).. Should try to make them more
        %consistent.
        
        
        
        %----------------------------------------------------------------
        % Switches for MODEL LWP, etc. data  
        %----------------------------------------------------------------
        
        %iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
        iadd_Conv_LWP=1;
        iadd_Conv_RWP=1;
        
        iuse_conv_for_LWP_to_TLWP_ratio=1; %whether to use the convective LWP and RWP in calculating the model LWP to LWP+RWP ratio.
        
        ilow_clear_only=0; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
        %mid or high cloud
        % N.B. - gets overruled and set to zero if irestrict_using_LWP_fraction=1    
        
        cloud_scene_selection = 'low_clear_only';
        cloud_scene_selection = 'MODIS COSP CTP screening';
        cloud_scene_selection = 'Model CTH screening';
        %cloud_scene_selection = 'none';
        
        model_CTH_thresh = 3.2e3; %CTH threshold in km
        
        CTP_thresh = 500e2; %CTP below which want to include clouds (Pa)
        CTP_thresh = 0; %CTP below which want to include clouds (Pa)
  
        irestrict_model_CF = 1;  %whether to restrict the model CFs to > 80%                                 
        
        imodis_cf = 1; %Whether to use the COSP MODIS liquid cloud fractions (low, mid and high), or the model ones for selecting low only scenes.
          % and also for screening for CF>80 etc. when using the CTH
          % method, or other methods for selecting low cloud scenes.

          
       % *** MAKE sure that also select the corresponding MODIS satellite
       % CF screening if doing model screening modis_CTH_screening)
        
        %----------------------------------------------------------------
        % Switches for SATELLITE LWP, etc. data
        %----------------------------------------------------------------
        
        modis_CTH_screening = '3.2km, CF>80';
        %modis_CTH_screening = '3.2km';
        %modis_CTH_screening = 'none';
        
        LWP_sat = 'AMSR-E';
        %LWP_sat = 'MODIS';
        
        %----------------------------------------------------------------
        
        %isave_plot_LWP_bias=0;
        fsize_latlon = 14;
        
    end
    
    
    if irestrict_using_LWP_fraction==1
        ilow_clear_only=0; %Use all points if using this for now
    end

    %% Set colourscale limits, etc. based on variable and switch choices.
    flag_str='';
    if iadd_RWP==1
        flag_str=[flag_str 'RWP'];
    end
    if iadd_Conv_LWP==1
        flag_str=[flag_str 'LWPconv'];
    end
    if iadd_Conv_RWP==1
        flag_str=[flag_str 'RWPconv'];
    end
    if irestrict_using_LWP_fraction==1
        flag_str=[flag_str 'restricted_LWPratio'];
    end
    

    switch plot_var_Sc
        case 'Ndat'
            clab01='Number of days';
            clab02='\DeltaN_{days}';
            clab01='Fraction of days';
            clab02='\Deltaf_{days}';
            caxis_range = [0 100];
            caxis_range_diff = [-100 100];
            caxis_range = [0 0.8];
            caxis_range_diff = [-0.2 0.2];
            var_Sc='Ndays_Sc_frac';
            flag_str='';
            
        case 'LWP'
            %defaults
            %var_Sc=['LWP' flag_str];
            var_Sc=['LWP'];
            clab_str='LWP';
            clab_units='(g m^{-2})';            
             switch bias_type
                case 'Absolute'
                    caxis_range_diff = [-100 100];
                    clab02='\DeltaLWP (g m^{-2})';
                case 'Percentage'
                    caxis_range_diff = [-100 100];
                    clab02='Bias (%)';
             end
             
            if iadd_RWP==1 
                clab_str=[clab_str '+RWP'];
            end
            
            if iadd_Conv_LWP==1 | iadd_Conv_RWP==1                
                if irestrict_using_LWP_fraction==1                    
                    caxis_range = [0 120];
                    switch bias_type
                        case 'Absolute'
                            caxis_range_diff = [-150 150];                            
                        case 'Percentage'
                            caxis_range_diff = [-150 150];                            
                    end
                else                    
                    caxis_range = [0 400];
                    switch bias_type
                        case 'Absolute'
                            caxis_range_diff = [-200 200];                            
                        case 'Percentage'
                            caxis_range_diff = [-200 200];                            
                    end
                end
                if iadd_Conv_LWP
                    clab_str=[clab_str '+LWP_{conv}'];
                end
                if iadd_Conv_RWP
                    clab_str=[clab_str '+RWP_{conv}'];
                end
            else
                if iadd_RWP==1
                    caxis_range = [0 230];
                    
                    
                    if irestrict_using_LWP_fraction==0                        
                        switch bias_type
                            case 'Absolute'
                                caxis_range_diff = [-50 50];
                            case 'Percentage'
                                caxis_range_diff = [-75 75];
                        end                   
                    end
                else
                    if irestrict_using_LWP_fraction==0
                        caxis_range = [0 150];
                        switch bias_type
                            case 'Absolute'
                                caxis_range_diff = [-50 50];
                            case 'Percentage'
                                caxis_range_diff = [-50 50];
                        end
                    else
                        caxis_range = [0 100];
                        switch bias_type
                            case 'Absolute'
                                caxis_range_diff = [-50 50];
                            case 'Percentage'
                                caxis_range_diff = [-150 150];
                        end
                    end
                end
                
            end
            
            clab01=[clab_str ' ' clab_units];
            
            
        case 'LWPic'
            clab01='LWP _{in-cloud} (g m^{-2})';
            
            %if iadd_Conv_LWP_RWP==1
            if iadd_Conv_LWP==1 | iadd_Conv_RWP==1                
                if irestrict_using_LWP_fraction==1
                    caxis_range = [0 350];
                else
                    caxis_range = [0 1000];
                end
            else
                if irestrict_using_LWP_fraction==1
                    caxis_range = [0 170];
                else
                    caxis_range = [0 350];
                end
            end
            switch bias_type
                case 'Absolute'
                    caxis_range_diff = [-100 100];
                    clab02='\DeltaLWP (g m^{-2})';
                case 'Percentage'
                    if irestrict_using_LWP_fraction==1
                        caxis_range_diff = [-100 100];
                    else
                        caxis_range_diff = [-50 50];
                    end
                    clab02='LWP_{in-cloud} bias (%)';
            end
            var_Sc='LWPic';
                        
            
        case 'Ratio of LWP to TLWP'
            caxis_range = [0 1];
            caxis_range = [0.6 1];
            caxis_range_diff = [-0.2 0.2];
            
            clab01='LWP to TLWP ratio';
            clab02='Bias';
            
            %bias_type = 'Percentage'; 
            %bias_type = 'Absolute'; 
            var_Sc='Ratio_LWP_TLWP';
            
        case 'COSP MODIS low CF estimate test'
            caxis_range = [0 1];
            caxis_range_diff = [-0.2 0.2];
            
            clab01='COSP MODIS liquid cloud fraction';
            clab02='Bias';
            
            bias_type = 'Absolute'; 
            var_Sc='COSP_MODIS_liqCF';
            
            flag_str='';
            
        case {'Nd'}
            caxis_range = [0 300];
            caxis_range_diff = [-100 100];
            
            clab01='Droplet Concentration (cm^{-3})';
            clab02='Bias (%)';
            
            bias_type = 'Percentage';
            var_Sc='Nd';
            
            LWP_sat = 'MODIS Nd';
            lwp_day_night = 'day';  
            
            flag_str='';
            
        case {'Nd COSP'}
            caxis_range = [0 600];
            caxis_range = [0 300];
            caxis_range_diff = [-200 200];
            caxis_range_diff = [-100 100];
            
            clab01='Droplet Concentration (cm^{-3})';
            clab02='Bias (%)';
            
            bias_type = 'Percentage';
            var_Sc='Nd_COSP';
            
            LWP_sat = 'MODIS Nd';
            lwp_day_night = 'day';    
            
            flag_str='';
            
    end
    
    %if irestrict_using_LWP_fraction==1
    %    var_Sc = [var_Sc 'restricted_LWPratio'];
    %end
    
    var_Sc = [var_Sc flag_str];

clear gca xlab Hc1s


    

    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=3; %no. rows, columns
    figure('position',scrsz);

   
    
%% ----  Model ---------------------------
    hs01=subplot(xsub,ysub,1);
    isub=isub+1; hs{isub}=hs01;
    gca = hs01;

   
        
%um_data 

switch cloud_scene_selection
    case 'low_clear_only'
        %if ilow_clear_only==1
        inds_clear_low_ALL = unique( [inds_PD_ALL{1}; inds_PD_ALL{2};] );
        
        % UM data
        % Do for all gridboxes, not just for box region as before
        if imodis_cf==1
            [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_modisCF_PD_ALL,mid_modisCF_PD_ALL,high_modisCF_PD_ALL,thresh_CF_states);
        else
            [forcing_vs_cloud_state_PD_ALL,std_forcing_vs_cloud_state_PD_ALL,freq_cloud_state_PD_ALL,Nii_cloud_state_PD_ALL,inds_PD_ALL,std_rel_norm_PD_ALL,cloud_states_PD_ALL] = calc_SW_forcing_cloud_fraction_combinations(forcing_ALL,low_CF_PD_ALL,mid_CF_PD_ALL,high_CF_PD_ALL,thresh_CF_states);
        end
        
        %Whereas using intersect here should make it so that e.g. PI needs to
        %be in the clear state AND PD in the clear state (but using the 4
        %combinations :- clear + clear, clear + low, etc.)
        % inds_CC = intersect( inds_PI{1}, inds_PD{1} );
        % inds_CL = intersect( inds_PI{1}, inds_PD{2} );
        % inds_LC = intersect( inds_PI{2}, inds_PD{1} );
        % inds_LL = intersect( inds_PI{2}, inds_PD{2} );
        
        %inds_PI_clear_low = unique( [inds_CC; inds_CL; inds_LC; inds_LL] );
        
    case 'MODIS COSP CTP screening'
        inds_clear_low_ALL = find(modis_cosp_CTP_PD >= CTP_thresh);
        
    case 'Model CTH screening'
        %inds_clear_low_ALL = find(max_cloud_height_PD <= model_CTH_thresh);
        %inds_clear_low_ALL = find(max_cloud_height_in_cloud_LWC_PD <= model_CTH_thresh);
        inds_clear_low_ALL = find(max_cloud_height_in_cloud_LWC_IWC_PD <= model_CTH_thresh);
        
        
    otherwise
        inds_clear_low_ALL = [1:prod(size(low_CF_PD_ALL))];
end

if imodis_cf==1
    %low_clear_CF_only = low_modisCF_PD_ALL(inds_clear_low_ALL);  %Low COSP
        %MODIS cloud fraction - not that reliable due to COSP putting
        %clouds too high in general. Probably need to add the mid-level
        %cloud.
    low_clear_CF_only = liqCF_COSP_MODIS_PD_ALL(inds_clear_low_ALL); %Liquid COSP MODIS CF - better to use this and select scenes
        %with low-only clouds (based on model CTH).
else
    low_clear_CF_only = low_CF_PD_ALL(inds_clear_low_ALL); 
end

i0=find(low_clear_CF_only>=-0.1);
i80=[1:prod(size(low_clear_CF_only))]; %default to all values - change below if set by the switches

if irestrict_model_CF==1
    i80=find(low_clear_CF_only>=0.8);
end


    
    LWP_LS_CONV = LWP_PD_ALL + Conv_LWP_PD_ALL; % Screening using the fraction of LWP to TLWP_ALL;
    TLWP_LS_CONV = LWP_LS_CONV + RWP_PD_ALL + Conv_RWP_PD_ALL;
    inan=find(TLWP_LS_CONV<0.1);
    TLWP_LS_CONV2 = TLWP_LS_CONV;
    %TLWP_LS_CONV2(inan) = 1; %Include these points as they are (almost) clear-sky
    %large-scale only values (not convective)
    TLWP_LS = LWP_PD_ALL + RWP_PD_ALL;
    inan_LS=find(TLWP_LS<0.1);
    if iuse_conv_for_LWP_to_TLWP_ratio==1
        LWP_fraction = LWP_LS_CONV ./ TLWP_LS_CONV2;
        LWP_fraction(inan)=1; %Include these points as they are (almost) clear-sky
    else
        %TLWP_LS (inan_LS) = 1; %Mistakenly had this set instead of LWP_fraction(inan_LS)=1 at one stage
        LWP_fraction = LWP_PD_ALL ./ TLWP_LS;
        LWP_fraction(inan_LS) = 1; %Include these points as they are (almost) clear-sky
    end
    %Need to calculate the above for later

    % Screening using the fraction of LWP to TLWP
if irestrict_using_LWP_fraction==1    
    i80 = find(LWP_fraction(inds_clear_low_ALL)>min_LWP_ratio); %i.e. mostly LWP    
end

 % Do time matching of model to satellite. 
     %Get the 3d array of local times for each location  in [lat,lon,time]
    %order.
    [LT] = local_times_from_UTC(HH_UM,gcm_Plon2D_UM);
    tol_hrs = 3; %tolerance for the matching (+/- N hours)
    switch lwp_day_night
        case 'average'
            inot_tol_LT=find(abs(LT-13.5)>tol_hrs & abs(LT-1.5)>tol_hrs);
        case 'day'
            inot_tol_LT=find(abs(LT-13.5)>tol_hrs);
        case 'night' 
            inot_tol_LT=find(abs(LT-1.5)>tol_hrs);
    end

add_str='';

switch plot_var_Sc   
    case {'LWP','LWPic'}
        dat_lwp = LWP_PD_ALL; % g/m2
                
        if iadd_RWP==1
            dat_lwp = dat_lwp + RWP_PD_ALL;
            add_str = ['including RWP'];
        else
            %subtitle_str=['Model'];
        end
        
        %if iadd_Conv_LWP_RWP==1
        if iadd_Conv_LWP==1
            dat_lwp = dat_lwp + Conv_LWP_PD_ALL;
            add_str = ['(' add_str ' & conv LWP)'];
        end
        if iadd_Conv_RWP==1
            dat_lwp = dat_lwp + Conv_RWP_PD_ALL;
            add_str = ['(' add_str ' & conv RWP)'];
        end
            %add_str = ['(' add_str ' & conv LWP+RWP)'];
        %end
        

        
        
        dat_modis_temp = NaN*ones(size(dat_lwp));
        dat_modis_temp(inds_clear_low_ALL(i80)) = dat_lwp(inds_clear_low_ALL(i80));
        
        
        
        
        
        dat_modis_temp(inot_tol_LT)=NaN; %discard nighttime values and ones outside of local time tol set above
        
        
        
        
        
        %Test using the data from the main code to check that they are similar
        [dat_modis, Ndat_modis, std_dev_modis] = meanNoNan(dat_modis_temp,3);
        
        switch plot_var_Sc
            case 'LWPic'
                cf = meanNoNan(low_calipsoCF_PD_ALL,3);                
                cf(cf<min_CF)=min_CF;
                dat_modis = dat_modis ./ cf;
                
                cf2 = low_calipsoCF_PD_ALL;
                cf2(cf2<min_CF)=min_CF;
                lwpic_model_all_times = dat_modis_temp ./ cf2;
                
        end
        
        %dat_modis = meanNoNan(low_clear_LWP_only_2d,3); %test
        
        add_str = [add_str ' day/night: ' lwp_day_night];
        
    case 'Ratio of LWP to TLWP'
        LWP_fraction_temp = NaN*ones(size(LWP_fraction));
        LWP_fraction_temp(inds_clear_low_ALL) = LWP_fraction(inds_clear_low_ALL); %don't use i80 restriction here since don't want to restrict to
        % fraction>0.99
        
        LWP_fraction_temp(inot_tol_LT)=NaN; %only use times when near 13 LT - will probably also want to do this comparison for the
        % case where we don't match with the MODIS data and so do for
        % both day and night (although would want to restrict to 01
        % and 13 LT).
        % Would do this in UM_ACSIS_global_SUBPLOT_commands_LWP.m
        
        [dat_modis, Ndat_modis, std_dev_modis] = meanNoNan(LWP_fraction_temp,3);
        
        
    case 'Ndat'
        
        %How to get the total available for the ratio?
        if iNtot_ratio_CF80_to_0 == 1
            dat_modis_temp2 = NaN*ones(size(dat_lwp));
            dat_modis_temp2(inds_clear_low_ALL(i0)) = 1;
        else
            dat_modis_temp2 = ones(size(dat_lwp));%use all data initially
        end
        
        dat_modis_temp2(inot_tol_LT)=NaN; %Don't include locations/times that are not in daylight for consistency with MODIS
        
        
        [temp, Ndat_tol, std_dev_tol] = meanNoNan(dat_modis_temp2,3);
        dat_modis=Ndat_modis./Ndat_tol;
        
    case 'Nd'
        dat_modis_temp = NaN*ones(size(dat_lwp));
        dat_modis_temp(inds_clear_low_ALL(i80)) = Nd_PD_ALL(inds_clear_low_ALL(i80))/1e6;  %convert to per cc
        dat_modis_temp(inot_tol_LT)=NaN; %discard nighttime values and ones outside of local time tol set above
        [dat_modis, Ndat_modis, std_dev_modis] = meanNoNan(dat_modis_temp,3);
        
        add_str = [add_str ' ' cloud_scene_selection ', irestrict model CF=' num2str(irestrict_model_CF) ', imodis cf=' num2str(imodis_cf)];
        
    case 'Nd COSP'
        dat_modis_temp = NaN*ones(size(dat_lwp));        
        %fad=1.0;
        %k=0.8;
        
        fad=0.7;
        k=0.88;
        
        CTT=273;
        CTT=263;
        

        [Nd_COSP,H,W,k,Q,cw]=MODIS_N_H_func(modis_cosp_tau_PD,modis_cosp_reff_PD,'calc',NaN,CTT,fad,k);
        dat_modis_temp(inds_clear_low_ALL(i80)) = Nd_COSP(inds_clear_low_ALL(i80));  %already in per cc
        dat_modis_temp(inot_tol_LT)=NaN; %discard nighttime values and ones outside of local time tol set above
        [dat_modis, Ndat_modis, std_dev_modis] = meanNoNan(dat_modis_temp,3);
        
        add_str = [add_str ' COSP ' cloud_scene_selection ', irestrict model CF=' num2str(irestrict_model_CF) ', imodis cf=' num2str(imodis_cf)];
        
        
    case 'COSP MODIS low CF estimate test'
        dat_modis_temp = NaN*ones(size(dat_lwp));
        dat_modis_temp(inds_clear_low_ALL(i80)) = liqCF_COSP_MODIS_PD_ALL(inds_clear_low_ALL(i80));
        
        [dat_modis, Ndat_modis, std_dev_modis] = meanNoNan(dat_modis_temp,3);
        
        %dat_modis = meanNoNan(liqCF_COSP_MODIS_PD_ALL,3);
        %dat_modis = meanNoNan(low_modisCF_PD_ALL + mid_modisCF_PD_ALL ,3);
end



    
    
    
    
    %Degrade UM data to that of CALIPSO (2x2 deg)
%    um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);    
    

    um_data = dat_modis; %save for later

    

    %Plotting Script
%    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %If plotting at UM native res
    %MODIS_ACSIS_LWP_global_vs_nest_quick_plot_commands_global    %
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)
    
    
%     ncols=55;
%     nlines=3;
%     istart_third = 2*ncols+1;
%     third_line = ['UM, mean=' num2str(Pmean,'%.2f')];
%     title_full = repmat(' ',[1 ncols*nlines]);
%     title_full(istart_third:istart_third+length(third_line)-1) = third_line;
%     title_full(1:length(subtitle_str)) = subtitle_str;
%     [tit_wrapped] = wrap_title_to_nlines(title_full,ncols,nlines);
    
    subtitle_str='Model';

    clear title_full       
    title_full{1}=[subtitle_str add_str ', mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str];
    title(title_full);
    
   
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
    %N.B. need to get a new handle each time want to change a different
    %colobar  - doesn't seem to allow multiple handles
    gca = hs01; Hc1s = find_peer_colorbars_of_an_axes(gca);    
    delete(Hc1s)  %just need one big colorbar for plots 1 and 2


%% ----  Sat data ---------------------------

    
    hs02=subplot(xsub,ysub,2);
    isub=isub+1; hs{isub}=hs02;
    gca = hs02;
    subtitle_str=LWP_sat;    
    
add_str ='';
    switch LWP_sat
        case 'AMSR-E'
            
            [date_str,date_num] = date_from_day_of_year_func(modis_loaded.daynum_timeseries3,modis_loaded.modisyear_timeseries3);
            
            %Time match the MODIS CTH data to the AMSRE days
                %Can just use when the CTH field is not NaN as a test of whether the CF>80
                %or CTH<3.2km criteria was met
                switch modis_CTH_screening
                    case '3.2km, CF>80'
                        [CTH_amsre] = time_match_data(amsre_matlab_time(time_inds_AMSRE),date_num,CTH.timeseries3);
                    case '3.2km';
                        [CTH_amsre] = time_match_data(amsre_matlab_time(time_inds_AMSRE),date_num,modis_loaded_CF0.CTH.timeseries3);
                    case 'none'
                        [CTH_amsre] = time_match_data(amsre_matlab_time(time_inds_AMSRE),date_num,modis_loaded_CF0_2.CTH.timeseries3);
                end
                                             
            %[CFliq_modis_amsre] = time_match_data(amsre_matlab_time(time_inds_AMSRE),date_num,Cloud_Fraction_Liquid.timeseries3);
            
            

            
            
            %qpcolor(amsre_LWP_Sc)
            %qpcolor(N_amsre_LWP_Sc)
            
            %Choose which AMSRE overpasses to use - do the same in the
            %model.
            switch lwp_day_night
                case 'average'
                    amsre_LWP_CTH = meanNoNan( 1e3*amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,1:2) ,4 ); %convert to g/m3 here
                    amsre_TLWP_CTH = meanNoNan( 1e3*amsre_tlwp(:,:,time_inds_AMSRE,1:2) ,4 ); %average over day and night overpasses
                    inan=[];
                case 'day'
                    amsre_LWP_CTH = 1e3*amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,1); %convert to g/m3 here
                    amsre_TLWP_CTH = 1e3*amsre_tlwp(:,:,time_inds_AMSRE,1);
                    [me,N_dat_Sc,stdev] = meanNoNan( CTH_amsre ,3);
                    
                    %Only do the MODIS-based screening when we are doing
                    %daytime AMSRE
                    
                    %When CTH data is NaN then the CTH<3.2km and CF>thresh_CF were
                    %not met
                    inan=find(isnan(CTH_amsre)==1);
                case 'night'
                    amsre_LWP_CTH = 1e3*amsre_dat.lwp_amsre(:,:,time_inds_AMSRE,2); %convert to g/m3 here
                    amsre_TLWP_CTH = 1e3*amsre_tlwp(:,:,time_inds_AMSRE,2);
                    inan=[];
            end
                      

                
                amsre_LWP_CTH = flipdim(amsre_LWP_CTH,1); %Flip the lat axis since AMSRE is oriented differently to MODIS
                amsre_LWP_CTH(inan)=NaN;
                
                
                amsre_TLWP_CTH = flipdim(amsre_TLWP_CTH,1); %Flip the lat axis since AMSRE is oriented differently to MODIS
                amsre_TLWP_CTH(inan)=NaN;
                
                lwp_ratio_amsre = amsre_LWP_CTH ./ amsre_TLWP_CTH;
                lwp_ratio_amsre(amsre_TLWP_CTH<min_LWP)=1; %set ratio to one when have near clear-sky, so that it is included.
                
                if irestrict_using_LWP_fraction==1                    
                    inan2 = find(lwp_ratio_amsre <= min_LWP_ratio | isnan(lwp_ratio_amsre)==1 );
                    amsre_LWP_CTH(inan2)=NaN;
                    amsre_TLWP_CTH(inan2)=NaN;                    
                end
                
                [amsre_LWP_Sc, N_amsre_LWP_Sc, stdev]  = meanNoNan( amsre_LWP_CTH ,3);                
                [amsre_TLWP_Sc, N_amsre_TLWP_Sc, stdev]  = meanNoNan( amsre_TLWP_CTH ,3);   

            
           %Degrade AMSR-E data to that of UM (2x2 deg for N96)            
                
                %flipdim requried here due to AMSRE data being flipped in
                %lat relative to its lat and lon arrays.
                switch plot_var_Sc
                    case 'Ndat'
                        
                        if iNtot_ratio_CF80_to_0==1
                            [me,Ntot_modis,stdev] = meanNoNan( CTH_amsre_CF0 ,3);
                        else
                            Ntot_modis = size(amsre_TLWP_CTH,3); %Assuming here that MODIS would be able to observe all times - however, not true due to some
                            %filtering for CTH, etc. So, trying ratio of CF>80 to
                            %CF>=0 for low altitude cloud only (mean CTH>3.2km).
                            %Think if there are no clouds then the CTH filtering
                            %does not apply? I.e., lets that gridpoint through?
                        end
                        
                        %dat_modis = griddata(flipdim(gcm_Plat2D_AMSRE_orig,1),flipdim(gcm_Plon2D_AMSRE_orig,1),N_amsre_TLWP_Sc./size(amsre_TLWP_CTH,3),gcm_Plat2D_UM,gcm_Plon2D_UM); %convert from kg/m2 to g/m2
                        dat_modis = griddata(flipdim(gcm_Plat2D_AMSRE_orig,1),flipdim(gcm_Plon2D_AMSRE_orig,1),N_dat_Sc./Ntot_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
                        %Better to use N_dat_Sc from CTH_amsre rather than
                        %N_amsre_TLWP_Sc since the latter depends on there
                        %also being an AMSRE LWP measurement
                        
                    case {'LWP','LWPic'}
                        if iadd_RWP==1 | iadd_Conv_RWP==1
                            dat_modis = griddata(flipdim(gcm_Plat2D_AMSRE_orig,1),flipdim(gcm_Plon2D_AMSRE_orig,1),amsre_TLWP_Sc,gcm_Plat2D_UM,gcm_Plon2D_UM);
                        else
                            dat_modis = griddata(flipdim(gcm_Plat2D_AMSRE_orig,1),flipdim(gcm_Plon2D_AMSRE_orig,1),amsre_LWP_Sc,gcm_Plat2D_UM,gcm_Plon2D_UM);
                        end
                        
                        switch plot_var_Sc
                            case 'LWPic'
                                cf_cal = meanNoNan(cllcalipso_monthly_AVERAGE(time_inds_CAL,:,:),1)/100;
                                %cf = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,cal_data,gcm_Plat2D_UM,gcm_Plon2D_UM);                                
                                cf = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,cf_cal,gcm_Plat2D_UM,gcm_Plon2D_UM);
                                cf(cf<min_CF)=min_CF;
                                dat_modis = dat_modis ./ cf;
                        end
                        
                
                        
                    case 'Ratio of LWP to TLWP'
                        lwp_ratio_amsre_mean  = meanNoNan( lwp_ratio_amsre ,3);
                        dat_modis = griddata(flipdim(gcm_Plat2D_AMSRE_orig,1),flipdim(gcm_Plon2D_AMSRE_orig,1),lwp_ratio_amsre_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);                        
                end
 
            
            %subtitle_str = [subtitle_str 'day/night setting:' lwp_day_night];
            
            add_str = [add_str ' day/night: ' lwp_day_night];
            
        case {'MODIS','MODIS Nd'}  
            switch LWP_sat
                case 'MODIS'
                    %Time match the MODIS CTH data to the AMSRE days
                    %Can just use when the CTH field is not NaN as a test of whether the CF>80
                    %or CTH<3.2km criteria was met
                    switch modis_CTH_screening
                        case '3.2km, CF>80'
                            %Degrade data to that of UM (2x2 deg for N96)
                            modis_LWP_time_period_mean = meanNoNan( modis_loaded.W_time3 .* modis_loaded.Cloud_Fraction_Liquid.timeseries3 ,3); %convert to grid-box mean
                            dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);                            
                        case '3.2km';
                            %Degrade data to that of UM (2x2 deg for N96)                            
                            modis_LWP_time_period_mean = meanNoNan( modis_loaded_CF0.W_time3 .* modis_loaded_CF0.Cloud_Fraction_Liquid.timeseries3 ,3); %convert to grid-box mean
                            dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
                            
                        case 'none'
                            %Degrade data to that of UM (2x2 deg for N96)                            
                            modis_LWP_time_period_mean = meanNoNan( modis_loaded_CF0_2.W_time3 .* modis_loaded_CF0_2.Cloud_Fraction_Liquid.timeseries3 ,3); %convert to grid-box mean
                            dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
                    end
                            
                            
                case 'MODIS Nd'
                    %Time match the MODIS CTH data to the AMSRE days
                    %Can just use when the CTH field is not NaN as a test of whether the CF>80
                    %or CTH<3.2km criteria was met
                    
                    if icorrect_MODIS_Nd_y2009==1                        
                        %Load the trend correction for y2009 vs y2000 - add to
                        %the MODIS data to make the obs more akin to 2009.
                        save_file_dNd_trend = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/dNd_trend_correction.mat'];
                        load(save_file_dNd_trend);                        
                    end
                    
                    switch modis_CTH_screening
                        case '3.2km, CF>80'
                            %Degrade data to that of UM (2x2 deg for N96)
                            modis_LWP_time_period_mean = meanNoNan( modis_loaded.Droplet_Number_Concentration_37.timeseries3 ,3);
                            dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);                            
                        case '3.2km';
                            %Degrade data to that of UM (2x2 deg for N96)
                            modis_LWP_time_period_mean = meanNoNan( modis_loaded_CF0.Droplet_Number_Concentration_37.timeseries3 ,3); 
                            dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);                            
                        case 'none'
                            %Degrade data to that of UM (2x2 deg for N96)
                            modis_LWP_time_period_mean = meanNoNan( modis_loaded_CF0_2.Droplet_Number_Concentration_37.timeseries3 ,3); 
                            dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_LWP_time_period_mean,gcm_Plat2D_UM,gcm_Plon2D_UM);
                    end
                    
                    if icorrect_MODIS_Nd_y2009==1
                        dat_modis = dat_modis + dNd_corr_UM_grid; trend_corr_str = ' corrected using linear trend to match y2000';
                        add_str = [add_str ' ' modis_CTH_screening trend_corr_str];
                    end
                    
                                                                                    
            end
                        
                    
    end
    
    
    sat_data = dat_modis; %back up data for future use perhaps


    %Run plotting script
%    MODIS_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)    
    %titlenam_driver = [subtitle_str];        
    %title([titlenam_driver ' mean=' num2str(Pmean,'%.2f')]);
    clear title_full
    title_full{1}=[subtitle_str ' ' add_str ', mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str];
    title(title_full);
    
    caxis(caxis_range);
    gca = hs02;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
    %N.B. need to get a new handle each time want to change a different
    %colobar  - doesn't seem to allow multiple handles
    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    

%     %make them all bigger - seem to have to increase the height rather than the
%     %width for this
%     new_height=0.25;
%     gca=hs01; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
%     gca=hs02; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
%     gca=hs03; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
% 
%     

  

%% ----  Difference plot ---------------------------

    
    hs03=subplot(xsub,ysub,3);
    isub=isub+1; hs{isub}=hs03;
    gca = hs03;
    subtitle_str='Bias UM minus satellite';

    abs_bias = um_data - sat_data; 
    prc_bias = 100*(um_data - sat_data) ./sat_data;
    prc_bias(sat_data<0.1)=NaN;
    switch bias_type
        case 'Absolute'
            dat_modis = abs_bias;
        case 'Percentage'
            dat_modis = prc_bias;          
    end
    

    
   
   if ibias_contour==1
       icontour_DRIVER=1;
       cont_dat = LWP_ratio_AMSRE;
       cont_ints=[0.9 0.9];
   end

   
    %Run plotting script
%    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)        
    %titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    %title([titlenam_driver ' mean=' num2str(Pmean,'%.2f')]);
    clear title_full
    title_full{1}=[subtitle_str , 'mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=['Model bias'];
    %N.B. - could then textwrap this to a certain number of columns
    title(title_full);
    
    caxis(caxis_range_diff);
    gca = hs03;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
    %N.B. need to get a new handle each time want to change a different
    %colobar  - doesn't seem to allow multiple handles
    Hc1s = find_peer_colorbars_of_an_axes(hs03);
   

    %adjust the position of the 2nd and 3rd plots to reduce space between the
    %subplots
    new_xpos=0.35;
    pos=get(hs02,'position');
    pos_orig=pos;
    pos(1)=new_xpos;
    set(hs02,'position',pos);
    
    new_xpos = 0.57;
    pos=get(hs03,'position');
    pos_orig=pos;
    pos(1)=new_xpos;
    set(hs03,'position',pos);


%Need to sort move the colourbars separately    
    


    %Make the 1st colorbar wider and move to below the plot
    Hc1s = find_peer_colorbars_of_an_axes(hs02);  
    pos=get(Hc1s,'position');
    pos(1)=0.12;
    pos(2)=0.24; cb_vert_pos = pos(2);
    pos(3)=0.4;
    set(Hc1s,'position',pos);
    
    
    %Adjust the 2nd colorbar
    Hc1s = find_peer_colorbars_of_an_axes(hs03);    
    pos=get(Hc1s,'position');
    pos(1)=0.5698;
    pos(2)=cb_vert_pos;
    pos(3)=0.17;
    set(Hc1s,'position',pos);

    %Change the colormap to blue to brown :-
    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1));    
%    dz_tit = 0.3;
    dz_tit = -0.7;    
    
    set(Hc1s,'XaxisLocation','bottom');
    xlab = xlabel(Hc1s,clab02,'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);    

    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    
    
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,[clab01],'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the label to move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  
    
    axes(hs01);
    pos_fig = get(gca,'position');


     if isave_plot_LWP_bias_Sc==1
%        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_LWP_bias_3panel_vs_satellite'];

        titlenam_driver=[var_Sc '_in_Sc_bias_3panel'];
        savename=[savedir_date titlenam_driver];
        
        % Run script to calculate biases, etc. for the different regions.
        iplot_box_on_map=0;  
        var_Latex=var_Sc;
        DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region
        
        %do the a,b,c, etc. labels
        label_subfigs
        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        
        
        

        

        
     end
    
    clear ioverride_LWP_Sc_choices 
catch error
    rethrow(error)
    clear ioverride_LWP_Sc_choices 
end