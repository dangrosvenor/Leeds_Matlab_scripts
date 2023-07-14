    clear gca xlab Hc1s hs
    
    isub=0; 
    
    aod_day_night = 'day';
    
    
%    isave_plot_CF = 1;
    
    fsize_latlon = 14;    
    caxis_range = [0 0.6];    

    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=3; %no. rows, columns
    figure('position',scrsz);
    
    bias_type = 'Absolute';
    bias_type = 'Percentage';
            
    

   
    
%% ----  Model ---------------------------
    hs01=subplot(xsub,ysub,1);
    isub=isub+1; hs{isub} = hs01;
    gca = hs01;
    
    % Do time matching of model to satellite. 
     %Get the 3d array of local times for each location  in [lat,lon,time]
    %order.
    [LT] = local_times_from_UTC(HH_UM,gcm_Plon2D_UM);
    tol_hrs = 3; %tolerance for the matching (+/- N hours)
    switch aod_day_night
        case 'average'
            inot_tol_LT=find(abs(LT-13.5)>tol_hrs & abs(LT-1.5)>tol_hrs);
        case 'day'
            inot_tol_LT=find(abs(LT-13.5)>tol_hrs);
        case 'night' 
            inot_tol_LT=find(abs(LT-1.5)>tol_hrs);
    end
    
    
    subtitle_str='Model';  
    dat_modis_temp = aod_total;
    dat_modis_temp(inot_tol_LT)=NaN; %discard nighttime values and ones outside of local time tol set above
    
    um_data = meanNoNan(dat_modis_temp,3);
    
    
    var_Latex='TotAOD'; %Used for both model and satellite.
    iplot_box_on_map=0;
    
           
    %Degrade UM data to that of CALIPSO (2x2 deg)
    %um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);            
    %dat_modis = um_cal_res;
    dat_modis = um_data;
    

    %Plotting Script
    var_UM = subtitle_str;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %If plotting at UM native res
%    CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global    %CALIPSO res
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error
    %title([titlenam_driver ' mean=' num2str(Pmean,'%.2g')]);    
    caxis(caxis_range);
    gca = hs01;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
    %N.B. need to get a new handle each time want to change a different
    %colobar  - doesn't seem to allow multiple handles
    gca = hs01; Hc1s = find_peer_colorbars_of_an_axes(gca);    
    delete(Hc1s)  %just need one big colorbar for plots 1 and 2


%% ----  MODIS sat data ---------------------------
    hs02=subplot(xsub,ysub,2);
    isub=isub+1; hs{isub} = hs02;
    gca = hs02;
    subtitle_str='MODIS';

    %cal_data = modis_aod_timeav;

          
    %Regrid to grid of UM
    %sat_data = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,cal_data,gcm_Plat2D_UM,gcm_Plon2D_UM);
    sat_data = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,modis_aod_timeav,gcm_Plat2D_UM,gcm_Plon2D_UM);                            
    dat_modis = sat_data;
    

    %Run plotting script
    %CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    var_UM = subtitle_str;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    %titlenam_driver = [subtitle_str];        
    %title([titlenam_driver ' mean=' num2str(Pmean,'%.2g')]);
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
    isub=isub+1; hs{isub} = hs03;
    gca = hs03;
    subtitle_str='Model bias';

%    sat_data = griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,cal_data,gcm_Plat2D_UM,gcm_Plon2D_UM);    
    
    %abs_bias = um_cal_res - cal_data;
    %prc_bias = 100 * abs_bias ./ cal_data;
    abs_bias = um_data - sat_data;
    prc_bias = 100 * abs_bias ./ sat_data;    
    bias_str='AOD bias';
    
    switch bias_type
        case 'Absolute'
            dat_modis = abs_bias;
            caxis_range_diff = [-0.3 0.3];              
        case 'Percentage'
            dat_modis = prc_bias;
            caxis_range_diff = [-100 100]; 
            bias_str = [bias_str ' (%)'];
    end
    
      
    
    
    
    

    %Run plotting script
    var_UM = subtitle_str;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    %CALIPSO_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    %titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    %title([titlenam_driver ' mean=' num2str(Pmean,'%.2g')]);
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
    xlab = xlabel(Hc1s,bias_str,'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);    

    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,['AOD'],'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  
    
    
if isave_plot_CF==1
        %savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_CF_bias_3panel_' cloud_alt];
        titlenam_driver=['global_AOD_bias_3panel'];
        savename=[savedir_date titlenam_driver];
        
        % Run script to calculate biases, etc. for the different regions.        
        col_str = 'k';
        axes(hs02); %switch to the satellite plot
        DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region        
        
        %do the a,b,c, etc. labels
        label_subfigs
        
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        %        close(gcf);
        
        
end 

