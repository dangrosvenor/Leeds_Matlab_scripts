clear gca xlab Hc1s hs

isave_plot_offline=1; %whether to save the plot at the 
fig_3panel_savename = 'trends_bar_plot_subfigs';
icoarse_grain=0;

igeneric_plot=1; %since the data has been coarse grained, so use GENERIC instead of UM

fsize_latlon = 14;  
%caxis_range = [-15 15]; 

cbar_title{1} = '%';
cbar_title{2} = '%';
cbar_title{3} = '%';
cbar_title{4} = '%';
cbar_title{5} = '%';
cbar_title{6} = '%';

var_Latex='SWest';

bias_type = 'Absolute';
bias_type = 'Percentage';

iplot_mgrid_lines_DRIVER = 1;
ioverride_ticks_DRIVER = 0; %keep the default tick locations
    



    subplotting=1; %added this functionality to MODIS_vert_pen_NdError_DRIVER_template_global_map.m
    xsub=2; ysub=3; %no. rows, columns
    figure('position',scrsz);
    isub=0;

   
    

   


%% ----  Sulphate change ---------------------------
    caxis_range=[-300 300];
    isub=isub+1;
    
    hs{isub} =subplot(xsub,ysub,isub);
    
    gca = hs{isub};
    subtitle_str='% Column Sulphate change';   
    
    
    dat_load = load([savedir_date 'Percentage_change_in_Column_H2SO4_mass_(PD_minus_PI).mat']);
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
    %title(title_full);
    title(subtitle_str);
    caxis(caxis_range);
    gca = hs{isub};
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width

 

%% ----  BC change ---------------------------
    caxis_range = [-300 300];      
    
    isub=isub+1;
    hs{isub}=subplot(xsub,ysub,isub);        
    gca = hs{isub}
    subtitle_str='% Column BC change';
    
    dat_load = load([savedir_date 'Percentage_change_in_Column_BC_mass_(PD_minus_PI).mat']);
    dat_modis = dat_load.plot_data;
             

    %Run plotting script
%    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)        
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    clear title_full       
    title_full{1}=[titlenam_driver ' mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str];
    %title(title_full);
    title(subtitle_str);
    caxis(caxis_range);
    gca = hs{isub}
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
    
%% ----  OM change ---------------------------
    caxis_range = [-100 100];    
    
    isub=isub+1;
    hs{isub}=subplot(xsub,ysub,isub);        
    gca = hs{isub};
    subtitle_str='% Column OM change';
    
    dat_load = load([savedir_date 'Percentage_change_in_Column_OM_mass_(PD_minus_PI).mat']);
    dat_modis = dat_load.plot_data;
             

    %Run plotting script
%    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)        
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    clear title_full       
    title_full{1}=[titlenam_driver ' mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str];
    %title(title_full);
    title(subtitle_str);
    caxis(caxis_range);
    gca = hs{isub};
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
%% ----  Dust (total) ---------------------------
    caxis_range = [-10 10];
    isub=isub+1;
    hs{isub}=subplot(xsub,ysub,isub);
    gca = hs{isub};    

    %um_data = 
    dat_load = load([savedir_date 'Percentage_change_in_Column_dust_mass_(PD_minus_PI).mat']);
    dat_modis = dat_load.plot_data;
    gcm_Plat2D_GENERIC = dat_load.lat_data;
    gcm_Plon2D_GENERIC = dat_load.lon_data;
    gcm_Plat2D_edges_GENERIC = dat_load.lat_data_edges;
    gcm_Plon2D_edges_GENERIC = dat_load.lon_data_edges;        
    daynum_timeseries3_GENERIC = daynum_timeseries3_UM;
    gcm_time_UTC_GENERIC = gcm_time_UTC_UM;
      
    subtitle_str=['% Column dust change'];

    
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
    title(subtitle_str);
    caxis(caxis_range);
    gca = hs{isub};
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
    
    
%% ----  sea-salt chagne ---------------------------
    caxis_range = [-7 7];     
    
    isub=isub+1;
    hs{isub}=subplot(xsub,ysub,isub);    
    gca = hs{isub};
    subtitle_str='% Column sea-salt change';
    
    dat_load = load([savedir_date 'Percentage_change_in_Column_sea-salt_mass_(PD_minus_PI).mat']);
    dat_modis = dat_load.plot_data;
             

    %Run plotting script
%    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)        
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    clear title_full       
    title_full{1}=[titlenam_driver ' mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str];
    %title(title_full);
    title(subtitle_str);
    caxis(caxis_range);
    gca = hs{isub};
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width   
    
%% ----  CCN chagne ---------------------------

iplot=0;
if iplot==1
    caxis_range = [-300 300];      
    
    hs06=subplot(xsub,ysub,6);
    isub=isub+1;
    hs{isub} = hs06;
    gca = hs06;
    subtitle_str='% CCN_{0.2%} change';
    
    dat_load = load([savedir_date 'Percentage_change_in_CCN_0.2pct_(PD_minus_PI).mat']);
    dat_modis = dat_load.plot_data;
             

    %Run plotting script
%    MODIS_ACSIS_Nd_global_vs_nest_quick_plot_commands_global
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global %main difference between scripts is which grid is used (e.g., UM, MODIS, CALIPSO, etc.)        
    titlenam_driver = [subtitle_str];    %So, we overestimate Nd due to this error    
    clear title_full       
    title_full{1}=[titlenam_driver ' mean=' num2str(Pmean,'%.2f')];
    title_full{2}='';
    title_full{3}=[subtitle_str];
    %title(title_full);
    title(subtitle_str);
    caxis(caxis_range);
    gca = hs06;
    increase_font_size_map_figures
    %move colorbar down a bit and expand to full width
end
    
    
%% Adjust plots and colourbar
%     %make them all bigger - seem to have to increase the height rather than the
%     %width for this
%     new_height=0.25;
%     gca=hs01; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
%     gca=hs02; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
%     gca=hs03; pos=get(gca,'position'); pos(4)=new_height; set(gca,'position',pos);
% 
%     




plot_width = 0.3;
dY_pos_ALL = -0.05;
for isub=1:5    
    pos=get(hs{isub},'position');    
    pos(3)=plot_width; pos(4)=plot_width;
    pos(2) = pos(2) + dY_pos_ALL;
    set(hs{isub},'position',pos);     
end

dX_pos = -0.065;
%dX_pos = -0.005;
for isub=1:5    
    pos=get(hs{isub},'position');
    pos(1) = pos(1) + dX_pos;    
    set(hs{isub},'position',pos);  
    
    %move cbars by dY_pos_ALL (will only adjust the x-pos for the 2nd and
    %3rd columns)
    Hc = find_peer_colorbars_of_an_axes(hs{isub});    
    pos = get(Hc,'position');    
    pos(2) = pos(2) + dY_pos_ALL;
    set(Hc,'position',pos);
end

%Now move the plots in first column relative to the first plot to reduce gaps
dX_pos = -0.06;
for isub=[2 5]
    pos=get(hs{isub},'position');
    pos(1) = pos(1) + dX_pos;    
    set(hs{isub},'position',pos);    
    
    %same for cbars
    Hc = find_peer_colorbars_of_an_axes(hs{isub});    
    pos = get(Hc,'position');
    pos(1) = pos(1) + dX_pos;    
    set(Hc,'position',pos);
end

%Now move the plots in first column relative to the first plot to reduce gaps
dX_pos = dX_pos*2;
for isub=[3]    
    pos=get(hs{isub},'position');
    pos(1) = pos(1) + dX_pos;    
    set(hs{isub},'position',pos);   
    
    %same for cbars
    Hc = find_peer_colorbars_of_an_axes(hs{isub});    
    pos = get(Hc,'position');
    pos(1) = pos(1) + dX_pos;    
    set(Hc,'position',pos);
end

   
%% sort the colourbars out
iadjust_colbars=0;
if iadjust_colbars==1
    
 %N.B. need to get a new handle each time want to change a different
    %colobar  - doesn't seem to allow multiple handles
    gca = hs03; Hc1s = find_peer_colorbars_of_an_axes(gca);    
    delete(Hc1s)  %just need one big colorbar for plots 1 and 2       
    
    
%Need to sort move the colourbars separately       

    %Make the 2nd colorbar wider and move to below the plot
    Hc1s = find_peer_colorbars_of_an_axes(hs02);    
    pos=get(Hc1s,'position');
    pos(1)=0.345;
    pos(2)=0.24; cb_vert_pos = pos(2);
    pos(3)=0.4;
    set(Hc1s,'position',pos);
    
    
    %Adjust the 2nd colorbar
    Hc1s = find_peer_colorbars_of_an_axes(hs01);    
    pos=get(Hc1s,'position');
    pos(1)=0.13;
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

    Hc1s = find_peer_colorbars_of_an_axes(hs02);
    
    set(Hc1s,'XaxisLocation','bottom');    
    xlab = xlabel(Hc1s,cbar_title{2},'fontsize',16);
    set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
    pos = get(xlab,'position');
    set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);  
    
    
end

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
        %DRIVER_calc_biases_for_regions %puts output in region_biases{i} for region_choice_multi{i} region
        
 
        
    end