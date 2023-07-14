clear gca xlab Hc1s hs

isave_plot_offline=1; %whether to save the plot at the end

fsize_latlon = 14;  
caxis_range = [0 100];

var_Latex='SWcalcs';

bias_type = 'Absolute';
bias_type = 'Percentage';

iplot_mgrid_lines_DRIVER = 1;
ioverride_ticks_DRIVER = 0; %keep the default tick locations
    



    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=3; %no. rows, columns
    figure('position',scrsz);
    isub=0;

   
    
%% ----  Model ---------------------------
    hs01=subplot(xsub,ysub,1);
    isub=isub+1;
    hs{isub} = hs01;
    gca = hs01;    

    %um_data = 
    dat_load = load([savedir_date 'SW_TOA_model_time-mean_calc.mat']);
    dat_modis = dat_load.plot_data;
    dat_modis(Land_mask_PD==1)=NaN;
      
    subtitle_str=['SWTOA estimated from model values'];

    
    %Degrade UM data to that of CALIPSO (2x2 deg)
%    um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);    
    

    um_data = dat_modis; %save for later

    

    %Plotting Script
%    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %If plotting at UM native res
    %MODIS_ACSIS_LWP_global_vs_nest_quick_plot_commands_global    %
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error
    clear title_full       
    title_full{1}=[titlenam_driver ' mean=' num2str(Pmean,'%.2f')];
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
    isub=isub+1;
    hs{isub} = hs02;
    gca = hs02;
    subtitle_str='SWTOA estimated from observed values';    
    
    dat_load = load([savedir_date 'SW_TOA_calculated_from_obs_mean_values_using_average_CF.mat']);
    dat_modis = dat_load.plot_data;


    subtitle_str = [subtitle_str];


    
    sat_data = dat_modis; %back up data for future use perhaps


    %Run plotting script
%    MODIS_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)    
    titlenam_driver = [subtitle_str];            
    clear title_full       
    title_full{1}=[titlenam_driver ' mean=' num2str(Pmean,'%.2f')];
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




    




    
    

%% ----  Bias plot ---------------------------
    caxis_range_diff = [-50 50];  
    
    hs03=subplot(xsub,ysub,3);
    isub=isub+1;
    hs{isub} = hs03;
    gca = hs03;
    subtitle_str='Reconstructed SWTOA bias';
    
    dat_load = load([savedir_date 'SW_TOA_perturbation_due_to_linear_combination_(pct).mat']);
    dat_modis = dat_load.plot_data;
             

    %Run plotting script
%    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)        
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    clear title_full       
    title_full{1}=[titlenam_driver ' mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str];
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
    xlab = xlabel(Hc1s,['Bias in ' var_name_units],'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);    

    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,'SW TOA (W m^{-2})','fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  
    
    %do the a,b,c, etc. labels
    label_subfigs
    
    

    
    
    if isave_plot_offline==1
        %        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_SW_bias_3panel_vs_CERES'];
        
        titlenam_driver=['SWTOA_offline_3panel'];
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        
        iplot_box_on_map=0;        
        % Run script to calculate biases, etc. for the different regions.
        %DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region
        
 
        
    end