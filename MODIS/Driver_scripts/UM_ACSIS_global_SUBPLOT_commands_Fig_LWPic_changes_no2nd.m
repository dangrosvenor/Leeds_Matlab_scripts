clear gca xlab Hc1s hs

isave_plot_offline=1; %whether to save the plot at the 
fig_3panel_savename = 'LWPic_changes_no2nd';
var_Latex='dLWPicNoSec';
icoarse_grain=1; M_coarse_grain=3; N_coarse_grain=3;
igeneric_plot=0;

fsize_latlon = 14;  
itwo_cbars=0;


switch fc_change_type
    case 'Absolute'
        cbar_title{1} = 'LWP_{ic} change (PD minus PI; g m^{-2})';
        caxis_range = [-100 100];
        caxis_range = [-50 50];
    case 'Percentage'
        cbar_title{1} = '% LWP_{ic} change';
        caxis_range = [-15 15];
        caxis_range = [-12.5 12.5];
end





iplot_mgrid_lines_DRIVER = 1;
ioverride_ticks_DRIVER = 0; %keep the default tick locations    

    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=1; ysub=2; %no. rows, columns
    figure('position',scrsz);
    isub=0;

   
    
%% ----  Standard runs with 2nd indirect effect ---------------------------
    hs01=subplot(xsub,ysub,1);
    isub=isub+1;
    hs{isub} = hs01;
    gca = hs01;    
    
    %dat_load = load([savedir_date 'Total_estimated_change_in_surface_SW.mat']);
    %dat_modis = dat_load.plot_dta;

    switch fc_change_type
        case 'Absolute'
            dat_modis = W1_mean - W0_mean; %var_UM = 'Absolute change in cloud fraction (PD minus PI)';
        case 'Percentage'
            dat_modis = 100*(W1_mean - W0_mean) ./ W0_mean; %var_UM = 'Percent change in cloud fraction, aerosol not affecting autoconversion (PD minus PI)';
    end

    gcm_Plat2D_GENERIC = dat_load.lat_data;
    gcm_Plon2D_GENERIC = dat_load.lon_data;
    gcm_Plat2D_edges_GENERIC = dat_load.lat_data_edges;
    gcm_Plon2D_edges_GENERIC = dat_load.lon_data_edges;        
    daynum_timeseries3_GENERIC = daynum_timeseries3_UM;
    gcm_time_UTC_GENERIC = gcm_time_UTC_UM;
    
%     icoarse_grain=1; M_coarse_grain=3; N_coarse_grain=3;
%     dat_modis = meanNoNan(direct_ALL,3);  
      
    subtitle_str=['Full model'];

    
    %Degrade UM data to that of CALIPSO (2x2 deg)
%    um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);    
    

% For stats will just calculate the % change between PI and PD for the no
% 2nd indirect run.
    %um_data = dat_modis; %save for later

    

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


%% ----  No 2nd indirect effect runs ---------------------------   
    hs02=subplot(xsub,ysub,2);
    isub=isub+1;
    hs{isub} = hs02;
    gca = hs02;
    subtitle_str='No aerosol effect on rain autconversion';    
    
    
    %dat_load = load([savedir_date 'Estimated_forcing_from_linear_sum_PI+PD_average.mat']);
    %dat_modis = dat_load.plot_data;
    switch fc_change_type
        case 'Absolute'
            dat_modis = W1_mean_no2nd - W0_mean_no2nd; %var_UM = 'Absolute change in cloud fraction (PD minus PI)';  
            %dat_modis = W1_mean_no2nd; %var_UM = 'Absolute change in cloud fraction (PD minus PI)';
        case 'Percentage'
            dat_modis = 100*((W1_mean_no2nd - W0_mean_no2nd) ./ W0_mean_no2nd); %var_UM = 'Percent change in cloud fraction, aerosol not affecting autoconversion (PD minus PI)';
    end
   
    gcm_Plat2D_GENERIC = dat_load.lat_data;
    gcm_Plon2D_GENERIC = dat_load.lon_data;
    gcm_Plat2D_edges_GENERIC = dat_load.lat_data_edges;
    gcm_Plon2D_edges_GENERIC = dat_load.lon_data_edges;        
    daynum_timeseries3_GENERIC = daynum_timeseries3_UM;
    gcm_time_UTC_GENERIC = gcm_time_UTC_UM;
    
%     icoarse_grain=1; M_coarse_grain=3; N_coarse_grain=3;
%     dat_modis = meanNoNan(indirect_ALL,3);  


    subtitle_str = [subtitle_str];

    %Save data for stats - want the % change from PI to PD for the no 2nd
    %indirect effect run
    um_data = W1_mean_no2nd;
    sat_data = W0_mean_no2nd; 


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
     
     
%% Move the subplots a bit    

    %adjust the position of the 2nd and 3rd plots to reduce space between the
    %subplots
    pos_orig=pos;
    new_xpos=0.45;
    pos=get(hs02,'position');    
    pos(1)=new_xpos;
    set(hs02,'position',pos);
    
%     new_xpos = 0.57;
%     pos=get(hs03,'position');
%     pos_orig=pos;
%     pos(1)=new_xpos;
%     set(hs03,'position',pos);


%% Sort the colourbars    
    
%Need to sort move the colourbars separately    
    
    %Make the 1st colorbar wider and move to below the plot
    Hc1s = find_peer_colorbars_of_an_axes(hs02);    
    pos=get(Hc1s,'position');
    pos(1)=0.12;
    pos(2)=0.16; cb_vert_pos = pos(2);
    pos(3)=0.61; %width of the 1st colourbar
    set(Hc1s,'position',pos);
    
    %Hc1s = find_peer_colorbars_of_an_axes(hs03); 
    
    if itwo_cbars==1
        %Adjust the 2nd colorbar
        
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
        xlab = xlabel(Hc1s,cbar_title{1},'fontsize',16);
        set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
        pos = get(xlab,'position');
        set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);
        
    else
        %delete(Hc1s);
    end

    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,cbar_title{1},'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  
    
    %do the a,b,c, etc. labels
    label_subfigs
    
    

    
    
    if isave_plot_offline==1
        %        savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_SW_bias_3panel_vs_CERES'];
        
        titlenam_driver=[fig_3panel_savename];
        savename=[savedir_date titlenam_driver];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        
        iplot_box_on_map=0;        
        % Run script to calculate biases, etc. for the different regions.
        DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region
        
 
        
    end
    
    icoarse_grain=0;
    
    